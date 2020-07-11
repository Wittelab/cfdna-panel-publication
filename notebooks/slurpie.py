from collections import OrderedDict
import pandas as pd
import glob
import seaborn as sns

from IPython.core.display import clear_output
from IPython.display import clear_output


## Accent Grey: #2e3436
## Muted Grey : #b2bec3
#flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71", "#fca50e"]
#flatui = ["#5deec5", "#85eceb", "#77bafc", "#a29dfb", "#fee9ab", "#fd9e4d", "#f8b1a2", "#fd7777", "#fb7ba8"]
flatui = ["#22cec8", "#1cb794", "#1a86e0", "#6c61e4", "#fcca75", "#fd9e4d", "#df7059", "#d43237", "#e64793"]
pal = sns.color_palette(flatui)


# filter_                     -- whether to apply the following filters after loading
# tolerable_alt_in_norm_frac  -- None to not tolerate any variants in the normal during filtering, 
#                                otherwise a normal variant frequency value between 0 and 1 for which 
#                                tumor variants present in the normal are tolerated
def read_gatk_tn(dir="TN/Study/wes-gatk/*.vcf", filter_=True, tolerable_alt_in_norm_frac=None):
    number_excluded = 0
    files = glob.glob(dir)
    foci = sorted(list(set(map(lambda k: k.split('/')[-1].split("_")[0], files))))
    fmap = OrderedDict(zip(foci, files))
    
    print files
    # Parse the data
    variants = None
    for n, (focus, f) in enumerate(fmap.iteritems(), 1):
        clear_output(wait=True)
        print "\rLoading %-16s (%d/%d) ..." % (focus, n, len(files))
        v = pd.read_csv(f, skiprows=57, header=None, sep='\t', error_bad_lines=False)
        print f

        v.insert(loc=0, column='Patient', value=focus.split('-')[0])
        v.insert(loc=1, column='Focus', value='-'.join(focus.split('-')[1:]).split('.')[0])
        variants = pd.concat([variants,v])
    
    variants.columns = [
                            "Patient",
                            "Focus",
                            "Chromosome", 
                            "Position", 
                            "ID", 
                            "Reference", 
                            "Variant", 
                            "Quality",
                            "Focus Filter", 
                            "Focus Information",
                            "Format",
                            "Tumor",
                            "Normal"
                       ]
    variants = variants.reset_index(drop=True)

    if filter_:
        number_excluded = variants.shape[0]
        if tolerable_alt_in_norm_frac!=None:
            aain = variants[variants['Focus Information'].apply(lambda k: 'alt_allele_in_normal'==k)]
            # Get the indices of variants that have a fraction lower than the acceptable limit 
            ok = (aain['Normal'].str.split(':')\
                 .apply(lambda k: float(k[1])/float(k[0])) \
                 <= tolerable_alt_in_norm_frac)\
                 .index
            # Set passing variants as OK
            variants['Focus Filter'].iloc[ok] = "PASS"
            print "Reaccepted # of variants: {0:,} (@<{1:.2f} in normal)".format(ok.shape[0], tolerable_alt_in_norm_frac)
        # Now filter for only passing variants 
        variants = variants[variants['Focus Filter']=="PASS"]
        # Calculate how many were removed
        number_excluded -= variants.shape[0]
        print "Filtered # of variants: {:,}".format(number_excluded)
        print "Resulting # of variants: {:,}".format(variants.shape[0])

    return variants



# filter_             -- whether to apply the following filters after loading
# tumor_thresh        -- BCBio threshold for the tumor (best alt likelihood)
# normal_thresh       -- BCBio threshold for normal (best ref likelihood)
# thresh_ratio        -- BCBio threshold for T/N frequency; avoids calling low frequency tumors also present at low frequency in normals
# generated_filters   -- a list of any in ['Focus Filter', 'LOD_PASS', 'FREQ_PASS'] for which to generate filters
def read_freebayes_tn(dir="TN/Study/wes-speed/normal-tumor/*.vcf", filter_=False, tumor_thresh=3.5, normal_thresh=3.5, thresh_ratio=2.7, generated_filters=['Focus Filter', 'LOD_PASS', 'FREQ_PASS']):
    if generated_filters is None and filter_ is True:
        print "If you'd like to filter, please specify the generated_filters list"
    number_excluded = 0
    files = glob.glob(dir)
    foci = sorted(list(set(map(lambda k: k.split('/')[-1].split("_")[0], files))))
    fmap = OrderedDict(zip(foci, files))
    
    # Parse the data
    variants = None
    for n, (focus, f) in enumerate(fmap.iteritems(), 1):
        clear_output(wait=True)
        print "\rLoading %-16s (%d/%d) ..." % (focus, n, len(files))
        v = pd.read_csv(f, skiprows=57, header=None, sep='\t', error_bad_lines=False)
        print f

        v.insert(loc=0, column='Patient', value=focus.split('-')[0])
        v.insert(loc=1, column='Focus', value='-'.join(focus.split('-')[1:]).split('.')[0])

        v.columns = [
                        "Patient",
                        "Focus",
                        "Chromosome", 
                        "Position", 
                        "ID", 
                        "Reference", 
                        "Variant", 
                        "Quality",
                        "Focus Filter", 
                        "Focus Information",
                        "Format",
                        "Normal",
                        "Tumor"
                   ]

        """
        ## LOD filter
        Call SOMATIC variants from tumor/normal calls, LOD_PASS.
        Assumes tumor/normal called with tumor first and normal second, as done in bcbio
        implementation.
        
        Uses MuTect like somatic filter based on implementation in speedseq:
        https://github.com/cc2qe/speedseq/blob/e6729aa2589eca4e3a946f398c1a2bdc15a7300d/bin/speedseq#L62
        Extracts the genotype likelihoods (GLs) from FreeBayes, which are like phred scores
        except not multiplied by 10.0 (https://en.wikipedia.org/wiki/Phred_quality_score).
        
        For tumors, we retrieve the best likelihood to not be reference (the first GL) and
        for normal, the best likelhood to be reference.
        
        After calculating the likelihoods, we compare these to thresholds to pass variants
        at tuned sensitivity/precision. Tuning done on DREAM synthetic 3 dataset evaluations.
        We also check that the frequency of the tumor exceeds the frequency of the normal by
        a threshold to avoid calls that are low frequency in both tumor and normal. This supports
        both FreeBayes and VarDict output frequencies.
        
        ## Frequency filter
        Ensure frequency of tumor to normal passes a reasonable threshold.
        Avoids calling low frequency tumors also present at low frequency in normals,
        which indicates a contamination or persistent error.
        Adds FREQ_PASS column. 
        """

        ## Check if a LOD filter is possible
        gl_index=None
        try:
            gl_index = v['Format'][0].split(":").index("GL")
        except ValueError:
            print "LOD check not possible (skipping...)"

        ## Check if a frequency filter is possible
        ao_index, ro_index = None, None  
        try:  
            ao_index = v['Format'][0].split(":").index("AO")
            ro_index = v['Format'][0].split(":").index("RO")
        except ValueError:
            print "Frequency check not possible (skipping...)"
            
            

        passes_lod = []
        tumor_lods = []
        normal_lods = []
        tn_freq = []
        passes_freq = []
        print "Evaluating filters..."
        for index, row in v.iterrows():
            ### LOD Filter
            if gl_index is not None and 'LOD_PASS' in generated_filters:
                try:
                    tumor_gls = [float(x) for x in row['Tumor'].split(":")[gl_index].split(",")]
                    tumor_lod = max(tumor_gls[i] - tumor_gls[0] for i in range(1, len(tumor_gls)))
                    #tumor_lod = min(tumor_gls[0] - tumor_gls[i] for i in range(1, len(tumor_gls)))
                # No GL information, no tumor call (so fail it)
                except IndexError:
                    tumor_lod = -1.0
                tumor_lods.append(tumor_lod)
                try:
                    normal_gls = [float(x) for x in row['Normal'].split(":")[gl_index].split(",")]
                    normal_lod = min(normal_gls[0] - normal_gls[i] for i in range(1, len(normal_gls)))
                    #normal_lod = min(normal_gls[i] - normal_gls[0] for i in range(1, len(normal_gls)))
                # No GL inofmration, no normal call (so pass it)
                except IndexError:
                    normal_lod = normal_thresh
                normal_lods.append(normal_lod)
                passes_lod.append(normal_lod >= normal_thresh and tumor_lod >= tumor_thresh)
            
            ### Frequency Filter
            if ao_index is not None and ro_index is not None and 'FREQ_PASS' in generated_filters:
                def _calc_freq(item):
                    try:
                        ao = sum([int(x) for x in item.split(":")[ao_index].split(",")])
                        ro = int(item.split(":")[ro_index])
                        freq = ao / float(ao + ro)
                    except (IndexError, ValueError, ZeroDivisionError):
                        freq = 0.0
                    return freq
                tumor_freq, normal_freq = _calc_freq(row['Tumor']), _calc_freq(row['Normal'])
                tn_freq.append((tumor_freq, normal_freq))
                passes_freq.append(normal_freq <= 0.001 or normal_freq <= tumor_freq / thresh_ratio)

        # Store results
        if 'LOD_PASS' in generated_filters:
            v['LOD_TUMOR'] = tumor_lods
            v['LOD_NORMAL'] = normal_lods
            v['LOD_PASS'] = passes_lod
        if 'FREQ_PASS' in generated_filters:
            v['TN_FREQS'] = tn_freq
            v['FREQ_PASS'] = passes_freq

        variants = pd.concat([variants,v])
    

    variants = variants.reset_index(drop=True)

    if filter_:
        number_excluded = variants.shape[0]
        if 'Focus Filter' in generated_filters:
            print "Filtering on PASS"
            variants = variants[variants['Focus Filter']=="PASS"]
        if 'LOD_PASS' in generated_filters:
            print "Filtering on LOD"
            variants = variants[variants['LOD_PASS']]
        if 'FREQ_PASS' in generated_filters:
            print "Filtering on Freq"
            variants = variants[variants['FREQ_PASS']]
        # Calculate how many were removed
        number_excluded -= variants.shape[0]
        print "Filtered # of variants: {:,}".format(number_excluded)
    print "Resulting # of variants: {:,}".format(variants.shape[0])

    return variants



# filter_             -- whether to apply the following filters after loading
# normal_thresh       -- BCBio threshold for normal (best ref likelihood)
# thresh_ratio        -- BCBio threshold for T/N frequency; avoids calling low frequency in normals
# generated_filters   -- a list of any in ['Focus Filter', 'LOD_PASS', 'FREQ_PASS'] for which to generate filters
def read_freebayes_mets(dir="TN/Study/mets-wes-speed/*.vcf", filter_=False, normal_thresh=3.5, thresh_ratio=2.7, generated_filters=['Focus Filter', 'LOD_PASS', 'FREQ_PASS']):
    if generated_filters is None and filter_ is True:
        print "If you'd like to filter, please specify the generated_filters list"
    number_excluded = 0
    files = glob.glob(dir)
    foci = sorted(list(set(map(lambda k: k.split('/')[-1].split("_")[0], files))))
    fmap = OrderedDict(zip(foci, files))
    
    # Parse the data
    variants = None
    for n, (focus, f) in enumerate(fmap.iteritems(), 1):
        clear_output(wait=True)
        print "\rLoading %-16s (%d/%d) ..." % (focus, n, len(files))
        v = pd.read_csv(f, skiprows=57, header=None, sep='\t', error_bad_lines=False)
        print f

        v.insert(loc=0, column='Patient', value=focus.split('-')[0])
        v.insert(loc=1, column='Focus', value='-'.join(focus.split('-')[1:]).split('.')[0])

        v.columns = [
                        "Patient",
                        "Focus",
                        "Chromosome", 
                        "Position", 
                        "ID", 
                        "Reference", 
                        "Variant", 
                        "Quality",
                        "Focus Filter", 
                        "Focus Information",
                        "Format",
                        "Normal",
                   ]
        """
        ## LOD filter
        Call SOMATIC variants from met normal calls, LOD_PASS.
        
        Uses MuTect like somatic filter based on implementation in speedseq:
        https://github.com/cc2qe/speedseq/blob/e6729aa2589eca4e3a946f398c1a2bdc15a7300d/bin/speedseq#L62
        Extracts the genotype likelihoods (GLs) from FreeBayes, which are like phred scores
        except not multiplied by 10.0 (https://en.wikipedia.org/wiki/Phred_quality_score).
        
        For met normals, the best likelhood to be reference is retreived
        
        After calculating the likelihoods, we compare these to thresholds to pass variants
        at tuned sensitivity/precision. Tuning done on DREAM synthetic 3 dataset evaluations.
        
        ## Frequency filter
        Ensure frequency variant is above a threhold
        Adds FREQ_PASS column. 
        """

        ## Check if a LOD filter is possible
        gl_index=None
        try:
            gl_index = v['Format'][0].split(":").index("GL")
        except ValueError:
            print "LOD check not possible (skipping...)"

        ## Check if a frequency filter is possible
        ao_index, ro_index = None, None  
        try:  
            ao_index = v['Format'][0].split(":").index("AO")
            ro_index = v['Format'][0].split(":").index("RO")
        except ValueError:
            print "Frequency check not possible (skipping...)"
            
            

        passes_lod = []
        normal_lods = []
        met_freq = []
        passes_freq = []
        print "Evaluating filters..."
        for index, row in v.iterrows():
            ### LOD Filter
            if gl_index is not None and 'LOD_PASS' in generated_filters:
                try:
                    normal_gls = [float(x) for x in row['Normal'].split(":")[gl_index].split(",")]
                    normal_lod = min(normal_gls[0] - normal_gls[i] for i in range(1, len(normal_gls)))
                    #normal_lod = min(normal_gls[i] - normal_gls[0] for i in range(1, len(normal_gls)))
                # No GL inofmration, no normal call (so pass it)
                except IndexError:
                    normal_lod = normal_thresh
                normal_lods.append(normal_lod)
                passes_lod.append(normal_lod >= normal_thresh)
            
            ### Frequency Filter
            if ao_index is not None and ro_index is not None and 'FREQ_PASS' in generated_filters:
                def _calc_freq(item):
                    try:
                        ao = sum([int(x) for x in item.split(":")[ao_index].split(",")])
                        ro = int(item.split(":")[ro_index])
                        freq = ao / float(ao + ro)
                    except (IndexError, ValueError, ZeroDivisionError):
                        freq = 0.0
                    return freq
                normal_freq = _calc_freq(row['Normal'])
                met_freq.append(normal_freq)
                passes_freq.append(normal_freq <= 0.001)

        # Store results
        if 'LOD_PASS' in generated_filters:
            v['LOD_NORMAL'] = normal_lods
            v['LOD_PASS'] = passes_lod
        if 'FREQ_PASS' in generated_filters:
            v['MET_FREQS'] = met_freq
            v['FREQ_PASS'] = passes_freq

        variants = pd.concat([variants,v])
    

    variants = variants.reset_index(drop=True)

    if filter_:
        number_excluded = variants.shape[0]
        if 'Focus Filter' in generated_filters:
            print "Filtering on PASS"
            variants = variants[variants['Focus Filter']=="PASS"]
        if 'LOD_PASS' in generated_filters:
            print "Filtering on LOD"
            variants = variants[variants['LOD_PASS']]
        if 'FREQ_PASS' in generated_filters:
            print "Filtering on Freq"
            variants = variants[variants['FREQ_PASS']]
        # Calculate how many were removed
        number_excluded -= variants.shape[0]
        print "Filtered # of variants: {:,}".format(number_excluded)
    print "Resulting # of variants: {:,}".format(variants.shape[0])

    return variants



# skiprows     -- number of rows to skip before the first variant
# filter_      -- whether to apply the following filters after loading
# pass_only    -- whether to exclude variants not meeting 'PASS' thresholds imposed by mageri (only when filter_ is True)
# min_freq     -- minimum variant frequency (0-100)
# max_freq     -- maximum variant frequency (0-100)
# min_coverage -- minimum total coverage (variant or not) of a position
# min_fam      -- the minimum number of families a variant should be present in
def read_mageri(dir="cfDNA/mageri/*.vcf*", skiprows=7042, filter_=True, pass_only=True, min_freq=0.00, max_freq=100, min_cov=2, min_fam=2):
    files = glob.glob(dir)
    patient = sorted(list(set(map(lambda k: k.split('/')[-1].split("-")[0], files))))
    fmap = OrderedDict(zip(patient, files))
    
    # Parse the data
    variants = None
    for n, (p, f) in enumerate(fmap.iteritems(), 1):
        clear_output(wait=True)
        print "\rLoading %-4s (%d/%d)" % (p, n, len(files))
        v = pd.read_csv(f, skiprows=skiprows, header=None, sep='\t', error_bad_lines=False)
        print f

        v.insert(loc=0, column='Patient', value=p)
        variants = pd.concat([variants,v])
    variants.columns = [
                            "Patient",
                            "Chromosome", 
                            "Position", 
                            "ID", 
                            "Reference", 
                            "Variant", 
                            "Quality",
                            "Filter", 
                            "Information",
                            "Format",
                            "cfDNA",
                    ]
    # Cleanup columns
    print "Harmonizing columns..."
    variants['Start'] = variants['Chromosome'].apply(lambda k: k.split(':')[1].split('-')[0])
    variants['End'] = variants['Chromosome'].apply(lambda k: k.split(':')[1].split('-')[1])
    variants['Position'] = variants['Start']
    variants['Chromosome'] = variants['Chromosome'].apply(lambda k: k.split(':')[0])
    variants['Total Coverage'] = variants['Information'].apply(lambda k: k.split(';')[0].strip("DP=")).astype(int)
    variants['Frequency'] = variants['Information'].apply(lambda k: k.split(';')[1].strip("AF=")).astype(float)*100
    variants['Non-Reference Coverage'] = (variants['Total Coverage'].astype(int) * variants['Frequency'].astype(float)/100).astype(int)

    if filter_:
        variants = filter_variants(variants, pass_only, min_freq, max_freq, min_cov, min_fam)
    print "Resulting # of variants: {:,}".format(variants.shape[0])
    return variants


# filter_      -- whether to apply the following filters after loading
# rare_only    -- Use only rare alleles as flagged by curio
# min_freq     -- minimum variant frequency (0-100)
# max_freq     -- maximum variant frequency (0-100)
# min_coverage -- minimum total family coverage (variant or not) of a position
# min_fam      -- minimum number of families in which a variant is present
def read_curio(dir="cfDNA/curio/loose/*.csv", filter_=True, rare_only=False, min_freq=0.00, max_freq=100, min_cov=2, min_fam=2):
    files = glob.glob(dir)
    patients = sorted(list(set(map(lambda k: k.split('/')[-1].split('-')[0], files))))
    fmap = OrderedDict(zip(patients, files))
    
    variants = None
    for n, (p,f) in enumerate(fmap.iteritems(), 1):
        clear_output(wait=True)
        print "\rLoading %-4s (%d/%d)" % (p, n, len(files))
        # Load Data
        fh = pd.read_csv(f, skiprows=28, sep=',', error_bad_lines=False)
        print f

        # Convert frequency to numeric 
        fh['Frequency'] = fh['Frequency'].str.strip('%').apply(float)
        fh['Filter'] = '.'
        fh.insert(loc=0, column='Patient', value=p)
        variants = pd.concat([variants, fh])
    
    if rare_only:
        variants = variants[variants['Type']=='Rare Allele']
    if filter_:
        variants = filter_variants(variants, False, min_freq, max_freq, min_cov, min_fam)
    print "Resulting # of variants: {:,}".format(variants.shape[0])
    return variants


# skiprows     -- number of rows to skip before the first variant
# filter_      -- whether to apply the following filters after loading
# pass_only    -- whether to exclude variants not meeting 'PASS' thresholds imposed by mageri (only when filter_ is True)
# min_freq     -- minimum variant frequency (0-100)
# max_freq     -- maximum variant frequency (0-100)
# min_coverage -- minimum total coverage (variant or not) of a position
# min_fam      -- the minimum number of families a variant should be present in
# min_qual     -- the minimum quality score that variants must meet to be included
def read_freebayes(dir="cfDNA/freebayes/*.vcf*", skiprows=153, filter_=True, pass_only=False, min_freq=0.00, max_freq=100, min_cov=2, min_fam=2, min_qual=20):
    number_excluded = 0
    files = glob.glob(dir)
    patient = sorted(list(set(map(lambda k: k.split('/')[-1].split("-")[0], files))))
    fmap = OrderedDict(zip(patient, files))
    
    # Parse the data
    variants = None
    for n, (p, f) in enumerate(fmap.iteritems(), 1):
        clear_output(wait=True)
        print "\rLoading %-4s (%d/%d)" % (p, n, len(files))
        v = pd.read_csv(f, skiprows=skiprows, header=None, sep='\t', error_bad_lines=False)
        print f,
        v.insert(loc=0, column='Patient', value=p)

        discovered = v.shape[0]
        print "Harmonizing..."
        v = v[~v[4].str.contains(',')]
        v['Non-Reference Coverage'] = v[7].apply(lambda k: k.split(';')[5].strip("AO=")).astype(int)
        v['Total Coverage'] = v[7].apply(lambda k: k.split(';')[7].strip("DP=")).astype(int)
        # Recalculate frequency (Don't trust AF from freebayes!)
        v['Frequency'] = (v['Non-Reference Coverage'].astype(float)/v['Total Coverage'])*100

        if filter_:
            v = v[v[5]>min_qual]
            v = filter_variants(v, pass_only, min_freq, max_freq, min_cov, min_fam, verbose=False)
        number_excluded += (discovered - v.shape[0])

        variants = pd.concat([variants,v])
    variants.columns = [
                            "Patient",
                            "Chromosome", 
                            "Position", 
                            "ID", 
                            "Reference", 
                            "Variant", 
                            "Quality", 
                            "Filter", 
                            "Information", 
                            "Format", 
                            "cfDNA", 
                            "Non-Reference Coverage", 
                            "Total Coverage", 
                            "Frequency", 
                    ]

    print "Resulting # of variants: {:,}".format(variants.shape[0])
    print "Filtered # of variants: {:,}".format(number_excluded)
    return variants


# skiprows     -- number of rows to skip before the first variant
# filter_      -- whether to apply the following filters after loading
# pass_only    -- whether to exclude variants not meeting 'PASS' thresholds imposed by mageri (only when filter_ is True)
# min_freq     -- minimum variant frequency (0-100)
# max_freq     -- maximum variant frequency (0-100)
# min_coverage -- minimum total coverage (variant or not) of a position
# min_fam      -- the minimum number of families a variant should be present in
# min_qual     -- the minimum quality score that variants must meet to be included
def read_vardict(dir="cfDNA/vardict/*.vcf*", skiprows=61, filter_=True, pass_only=True, min_freq=0.00, max_freq=100, min_cov=2, min_fam=2, min_qual=20):
    files = glob.glob(dir)
    patient = sorted(list(set(map(lambda k: k.split('/')[-1].split("-")[0], files))))
    fmap = OrderedDict(zip(patient, files))
    
    # Parse the data
    variants = None
    for n, (p, f) in enumerate(fmap.iteritems(), 1):
        clear_output(wait=True)
        print "\rLoading %-4s (%d/%d)" % (p, n, len(files))
        v = pd.read_csv(f, skiprows=skiprows, header=None, sep='\t', error_bad_lines=False)
        print f

        v.insert(loc=0, column='Patient', value=p)
        variants = pd.concat([variants,v])
    variants.columns = [
                            "Patient",
                            "Chromosome", 
                            "Position", 
                            "ID", 
                            "Reference", 
                            "Variant", 
                            "Quality",
                            "Filter", 
                            "Information",
                            "Format",
                            "cfDNA",
                    ]
    # Cleanup columns
    print "Harmonizing columns..."
    variants['Total Coverage'] = variants['Information'].apply(lambda k: k.split(';')[2].strip("DP=")).astype(int)
    variants['Frequency'] = variants['Information'].apply(lambda k: k.split(';')[4].strip("AF=")).astype(float)*100
    variants['Non-Reference Coverage'] = variants['Information'].apply(lambda k: k.split(';')[3].strip("VD=")).astype(int)

    if filter_:
        variants = variants[variants['Quality']>min_qual]
        variants = filter_variants(variants, pass_only, min_freq, max_freq, min_cov, min_fam)
    print "Resulting # of variants: {:,}".format(variants.shape[0])
    return variants


def read_coverages(dir="cfDNA/coverages/*.csv"):
    files = glob.glob(dir)
    patients = sorted(list(set(map(lambda k: k.split('/')[-1].split('-')[0], files))))
    fmap = OrderedDict(zip(patients, files))
    
    coverages = None
    for n, (p,f) in enumerate(fmap.iteritems(), 1):
        clear_output(wait=True)
        print "\rLoading %-4s (%d/%d) ..." % (p, n, len(files))
        # Load Data
        fh = pd.read_csv(f, skiprows=35, sep=',', error_bad_lines=False)
        print f

        fh.insert(loc=0, column='Patient', value=p)
        coverages = pd.concat([coverages, fh])
    coverages['Coverage'] = coverages['Average Coverage Depth (Families)'] * \
    coverages['UMI/UMT Avg Family Size (Reads)']

    # Parse the actual positions
    coverages['Chromosome'] = coverages['Feature Location'].apply(lambda k: k.split(':')[0])
    coverages['Start'] = coverages['Feature Location'].apply(lambda k: k.split(':')[1].split('-')[0]).astype(int)
    coverages['End'] = coverages['Feature Location'].apply(lambda k: k.split(':')[1].split('-')[1]).astype(int)
    return coverages


def read_num_reads(file="cfDNA/sample_reads.tsv"):
    return pd.read_csv(file, sep="\t", header=None, index_col=0, names=["Reads"]).rename_axis('Patient')


# pass_only    -- only return variants that have 'PASS' for the filter
# min_freq     -- minimum variant frequency (0-100)
# max_freq     -- maximum variant frequency (0-100)
# min_coverage -- minimum total coverage (variant or not) of a position
# min_fam      -- the minimum number of families a variant should be present in
def filter_variants(variants, pass_only=False, min_freq=0.00, max_freq = 100, min_cov=2, min_fam=2, verbose=True):
    number_excluded = variants.shape[0]
    # Keep only those variants that have a PASS for filter
    if pass_only:
        variants = variants[variants['Filter']=='PASS']
    # Remove variants with frequency below min_freq threshold
    variants = variants[variants['Frequency']>=min_freq]
    # Remove variants with frequency above max_freq threshold
    variants = variants[variants['Frequency']<=max_freq]
    # Remove variants with low total read/family coverage
    variants = variants[variants['Total Coverage']>=min_cov]
    # Remove variants with few supporting variant reads/families
    variants = variants[variants['Non-Reference Coverage']>=min_fam]
    number_excluded -= variants.shape[0]
    if verbose:
        print "Filtered # of variants: {:,}".format(number_excluded)
    return variants 


# variants should have a 'Patient', 'Chromosome', and 'Position' column
# coverags should have a 'Patient', 'Chromosome', 'Start', 'End' and 'Coverage' column
# Will return the variants for which there is an intersection with a coverage region, and the coverage value
def get_variant_coverage(variants, coverages):
    found = None
    for patient in set(variants['Patient']):
        these_ = variants[variants['Patient']==patient]
        those_ = coverages[coverages['Patient']==patient]
        for chr_ in set(these_['Chromosome']):
            print "Processing... %s, %s" % (patient, chr_)
            clear_output(wait=True)
            this_ = these_[these_['Chromosome']==chr_]
            that_ = those_[those_['Chromosome']==chr_][['Start','End','Coverage']]
            fs = this_.merge(that_, how='inner', left_on='Position', right_on='Start')
            fe = this_.merge(that_, how='inner', left_on='Position', right_on='End')
            found = pd.concat([found, fs, fe])
    found = found.drop_duplicates()
    return found