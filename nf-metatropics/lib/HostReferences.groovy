import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.util.zip.GZIPInputStream
import groovy.transform.Immutable

/**
 * Resolve host background references from a friendly keyword, optionally downloading into <repo>/Databases.
 *
 * This is intended to keep the existing behaviour:
 * - If users set params.Human_host_fasta / params.Other_host_fasta explicitly → no changes.
 * - If users set params.Host (e.g. "human") and the matching *_host_fasta param is unset →
 *   download/cache the FASTA into Databases/ and set the correct param path.
 */
class HostReferences {

    @Immutable
    static class HostSpec {
        String key
        String group              // "human" → Human_host_fasta, everything else → Other_host_fasta
        String relativeDir        // under <repo>/Databases
        String fastaBasename      // file name after download (uncompressed)
        List<String> urls         // mirror list, try in order
        String sha256             // optional integrity check for downloaded bytes (compressed or uncompressed, see implementation)
    }

    // Central registry for keyword hosts.
    // Add more entries here as you grow the set of supported hosts.
    static final Map<String, HostSpec> REGISTRY = [
        // Human: CHM13 v2.0 (T2T)
        // Mirrors: primary is the T2T repo release FASTA (gz).
        // (If you prefer different sources, swap the URL list.)
        'human': new HostSpec(
            'human',
            'human',
            'Human',
            'chm13v2.0.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz'
            ],
            null
        ),

        // Pan troglodytes: NCBI RefSeq assembly (gz).
        // Note: exact assembly choice is a policy decision; this is a sensible default.
        'pan': new HostSpec(
            'pan',
            'other',
            'Pan',
            'pan_troglodytes.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz'
            ],
            null
        ),

        // Aedes aegypti (mosquito)
        'aedes': new HostSpec(
            'aedes',
            'other',
            'Aedes',
            'aedes_aegypti.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna.gz'
            ],
            null
        ),

        // Anopheles gambiae (mosquito)
        'anopheles': new HostSpec(
            'anopheles',
            'other',
            'Anopheles',
            'anopheles_gambiae.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/575/GCA_000005575.1_AgamP3/GCA_000005575.1_AgamP3_genomic.fna.gz'
            ],
            null
        ),

        // Culex quinquefasciatus (mosquito)
        'culex': new HostSpec(
            'culex',
            'other',
            'Culex',
            'culex_quinquefasciatus.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/732/765/GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal/GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna.gz'
            ],
            null
        ),

        // Macaca mulatta (rhesus macaque)
        'macaque': new HostSpec(
            'macaque',
            'other',
            'Macaque',
            'macaca_mulatta.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.fna.gz'
            ],
            null
        ),

        // Gorilla gorilla gorilla (western lowland gorilla)
        'gorilla': new HostSpec(
            'gorilla',
            'other',
            'Gorilla',
            'gorilla_gorilla.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/281/585/GCF_029281585.1_NHGRI_mGorGor1-v1.1-0.2.freeze_pri/GCF_029281585.1_NHGRI_mGorGor1-v1.1-0.2.freeze_pri_genomic.fna.gz'
            ],
            null
        ),

        // Pongo abelii (Sumatran orangutan)
        'orangutan': new HostSpec(
            'orangutan',
            'other',
            'Orangutan',
            'pongo_abelii.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/655/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna.gz'
            ],
            null
        ),

        // Sus scrofa (pig)
        'pig': new HostSpec(
            'pig',
            'other',
            'Pig',
            'sus_scrofa.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Sus_scrofa/latest_assembly_versions/GCF_000003025.6_Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic.fna.gz'
            ],
            null
        ),

        // Bos taurus (cow)
        'cow': new HostSpec(
            'cow',
            'other',
            'Cow',
            'bos_taurus.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz'
            ],
            null
        ),

        // Mus musculus (mouse)
        'mouse': new HostSpec(
            'mouse',
            'other',
            'Mouse',
            'mus_musculus.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz'
            ],
            null
        ),

        // Canis lupus familiaris (dog)
        'dog': new HostSpec(
            'dog',
            'other',
            'Dog',
            'canis_lupus_familiaris.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Canis_lupus_familiaris/latest_assembly_versions/GCF_011100685.1_UU_Cfam_GSD_1.0/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.fna.gz'
            ],
            null
        ),

        // Felis catus (cat)
        'cat': new HostSpec(
            'cat',
            'other',
            'Cat',
            'felis_catus.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Felis_catus/latest_assembly_versions/GCF_018350175.1_F.catus_Fca126_mat1.0/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna.gz'
            ],
            null
        ),

        // Rattus norvegicus (rat)
        'rat': new HostSpec(
            'rat',
            'other',
            'Rat',
            'rattus_norvegicus.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_genomic.fna.gz'
            ],
            null
        ),

        // Gallus gallus (chicken)
        'chicken': new HostSpec(
            'chicken',
            'other',
            'Chicken',
            'gallus_gallus.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz'
            ],
            null
        ),

        // Camelus dromedarius (dromedary camel)
        'camel': new HostSpec(
            'camel',
            'other',
            'Camel',
            'camelus_dromedarius.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/803/125/GCF_000803125.2_CamDro3/GCF_000803125.2_CamDro3_genomic.fna.gz'
            ],
            null
        ),

        // Capra hircus (goat)
        'goat': new HostSpec(
            'goat',
            'other',
            'Goat',
            'capra_hircus.refseq.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/704/415/GCF_001704415.2_ARS1.2/GCF_001704415.2_ARS1.2_genomic.fna.gz'
            ],
            null
        ),

        // Rhinolophus ferrumequinum (greater horseshoe bat) — representative bat host
        'bat': new HostSpec(
            'bat',
            'other',
            'Bat',
            'rhinolophus_ferrumequinum.genbank.fa',
            [
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/108/255/GCA_014108255.1_mRhiFer1.p/GCA_014108255.1_mRhiFer1.p_genomic.fna.gz'
            ],
            null
        )
    ].asImmutable()

    static Map resolve(params, log, baseDir) {
        def requested = normaliseHostParam(params.Host)
        if (!requested) {
            return [ human: null, other: null ]
        }

        // Don't download large references when printing help/version
        if (params.help || params.version) {
            return [ human: null, other: null ]
        }

        // If user explicitly provided paths, preserve them; Host only fills gaps.
        // Also respect legacy aliases (fasta / host_fasta), which main.nf already promotes reliably.
        def needsHuman = !params.Human_host_fasta && !params.fasta && requested.contains('human')
        def otherHosts = requested.findAll { it != 'human' }
        def needsOther = !params.Other_host_fasta && !params.host_fasta && otherHosts.size() > 0

        if (!needsHuman && !needsOther) {
            return [ human: null, other: null ]
        }

        Path dbRoot = Path.of(baseDir.toString(), 'Databases')
        Files.createDirectories(dbRoot)

        Path humanPath = null
        Path otherPath = null

        if (needsHuman) {
            def spec = REGISTRY['human']
            humanPath = ensureHostFasta(dbRoot, spec, log)
            log.info "Resolved Host='human' to Human_host_fasta: ${humanPath}"
        }

        if (needsOther) {
            def resolvedSpecs = otherHosts.collect { key ->
                def spec = REGISTRY[key]
                if (!spec) {
                    def supported = REGISTRY.keySet().sort().join(', ')
                    throw new IllegalArgumentException("Unsupported Host value '${key}'. Supported: ${supported}. " +
                            "Alternatively, set --Other_host_fasta / --Human_host_fasta explicitly.")
                }
                return spec
            }

            def resolvedPairs = resolvedSpecs.collect { spec -> [ spec.key, ensureHostFasta(dbRoot, spec, log) ] }
            resolvedPairs = resolvedPairs.unique { it[0] }.sort { a, b -> a[0] <=> b[0] }
            def resolvedKeys = resolvedPairs.collect { it[0] }
            def resolvedFastas = resolvedPairs.collect { it[1] }

            Path otherFasta
            if (resolvedFastas.size() == 1) {
                otherFasta = resolvedFastas[0]
            } else {
                otherFasta = ensureMergedOtherHostFasta(dbRoot, resolvedKeys, resolvedFastas)
                log.info "Merged non-human hosts (${resolvedKeys.join(', ')}) into Other_host_fasta: ${otherFasta}"
            }

            otherPath = otherFasta
            log.info "Resolved Host (${resolvedKeys.join(', ')}) to Other_host_fasta: ${otherPath}"
        }

        return [ human: humanPath?.toString(), other: otherPath?.toString() ]
    }

    private static List<String> normaliseHostParam(hostParam) {
        if (!hostParam) return []
        def items = []
        if (hostParam instanceof List) {
            items = hostParam
        } else {
            items = hostParam.toString().split(/[,\s]+/).findAll { it }
        }
        return items.collect { it.toString().trim().toLowerCase() }.findAll { it }.unique()
    }

    private static Path ensureHostFasta(Path dbRoot, HostSpec spec, log) {
        Path hostDir = dbRoot.resolve(spec.relativeDir)
        Files.createDirectories(hostDir)
        Path targetFasta = hostDir.resolve(spec.fastaBasename)
        if (Files.exists(targetFasta) && Files.size(targetFasta) > 0) {
            return targetFasta
        }

        // Download into a temp file first, then move into place.
        Path tmp = hostDir.resolve("${spec.fastaBasename}.download.tmp")
        Files.deleteIfExists(tmp)

        def lastErr = null
        for (String url : spec.urls) {
            try {
                log.info "Downloading host FASTA for '${spec.key}' from ${url}"
                downloadUrl(url, tmp)
                // If gzipped, decompress into target; otherwise move.
                if (isGzip(tmp)) {
                    gunzip(tmp, targetFasta)
                    Files.deleteIfExists(tmp)
                } else {
                    Files.move(tmp, targetFasta, StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE)
                }
                if (Files.size(targetFasta) == 0) {
                    throw new IOException("Downloaded FASTA is empty at ${targetFasta}")
                }
                return targetFasta
            } catch (Exception e) {
                lastErr = e
                Files.deleteIfExists(tmp)
                Files.deleteIfExists(targetFasta)
                log.warn "Failed downloading '${spec.key}' from ${url}: ${e.message}"
            }
        }
        throw new RuntimeException("Failed to download host FASTA for '${spec.key}'. Tried: ${spec.urls.join(', ')}", lastErr)
    }

    private static Path ensureMergedOtherHostFasta(Path dbRoot, List<String> keys, List<Path> fastaPaths) {
        Path mergedDir = dbRoot.resolve('_merged')
        Files.createDirectories(mergedDir)

        def safeKeys = keys.collect { it.toString().trim().toLowerCase() }.findAll { it }.unique().sort()
        def basename = "OtherHost__${safeKeys.join('+')}.fa"
        Path merged = mergedDir.resolve(basename)
        if (Files.exists(merged) && Files.size(merged) > 0) {
            return merged
        }

        Path tmp = mergedDir.resolve("${basename}.tmp")
        Files.deleteIfExists(tmp)

        tmp.toFile().withOutputStream { os ->
            fastaPaths.eachWithIndex { Path p, int idx ->
                if (idx > 0) os.write('\n'.getBytes('UTF-8'))
                p.toFile().withInputStream { is -> is.transferTo(os) }
            }
        }

        if (Files.size(tmp) == 0) {
            throw new IOException("Merged host FASTA is empty at ${tmp}")
        }
        Files.move(tmp, merged, StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE)
        return merged
    }

    private static void downloadUrl(String urlStr, Path out) {
        URL url = new URL(urlStr)
        url.openConnection().with { conn ->
            conn.setConnectTimeout(30_000)
            conn.setReadTimeout(300_000)
            conn.setRequestProperty('User-Agent', 'metatropics-nextflow')
            conn.connect()
            out.toFile().withOutputStream { os ->
                conn.getInputStream().withStream { is ->
                    is.transferTo(os)
                }
            }
        }
    }

    private static boolean isGzip(Path path) {
        path.toFile().withInputStream { is ->
            byte[] magic = new byte[2]
            if (is.read(magic) != 2) return false
            return (magic[0] == (byte)0x1f) && (magic[1] == (byte)0x8b)
        }
    }

    private static void gunzip(Path gzPath, Path outPath) {
        gzPath.toFile().withInputStream { fis ->
            new GZIPInputStream(fis).withStream { gis ->
                outPath.toFile().withOutputStream { os ->
                    gis.transferTo(os)
                }
            }
        }
    }
}

