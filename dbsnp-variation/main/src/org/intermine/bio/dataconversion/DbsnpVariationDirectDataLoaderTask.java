package org.intermine.bio.dataconversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;

import org.intermine.model.InterMineObject;
import org.intermine.model.bio.BioEntity;
import org.intermine.model.bio.SequenceFeature;
import org.intermine.model.bio.Chromosome;
import org.intermine.model.bio.SequenceAlteration;
import org.intermine.model.bio.AliasName;
import org.intermine.model.bio.Consequence;
import org.intermine.model.bio.ConsequenceType;
import org.intermine.model.bio.DataSet;
import org.intermine.model.bio.DataSource;
import org.intermine.model.bio.Gene;
import org.intermine.model.bio.Location;
import org.intermine.model.bio.Transcript;
import org.intermine.model.bio.Ontology;
import org.intermine.model.bio.Organism;
import org.intermine.model.bio.SNP;
import org.intermine.model.bio.Indel;
import org.intermine.model.bio.TandemRepeat;
import org.intermine.model.bio.Substitution;
import org.intermine.model.bio.SOTerm;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.objectstore.proxy.ProxyReference;
import org.intermine.objectstore.query.ConstraintSet;
import org.intermine.metadata.ConstraintOp;
import org.intermine.metadata.StringUtil;
import org.intermine.objectstore.query.ContainsConstraint;
import org.intermine.objectstore.query.QueryClass;
import org.intermine.objectstore.query.QueryField;
import org.intermine.objectstore.query.QueryObjectReference;
import org.intermine.objectstore.query.Query;
import org.intermine.objectstore.query.QueryValue;
import org.intermine.objectstore.query.Results;
import org.intermine.objectstore.query.ResultsRow;
import org.intermine.objectstore.query.SimpleConstraint;
import org.intermine.task.FileDirectDataLoaderTask;
import org.intermine.util.FormattedTextParser;

public class DbsnpVariationDirectDataLoaderTask extends FileDirectDataLoaderTask {

    private static final String DATASET_TITLE = "dbSNP (build 147) variants and variant annotations";
    private static final String DATA_SOURCE_NAME = "ensembl dbSNP";
    private static final Logger LOG = Logger.getLogger(DbsnpVariationDirectDataLoaderTask.class);
    private Integer taxonId = 9913;
    final static String[] expectedHeaders = {
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO"
    };
    private static final String VARIANT_ANNOTATION_SOURCE = "Ensembl VEP";
    private static final ArrayList < String > FUNCTION_CLASS_TO_IGNORE = new ArrayList < String > (Arrays.asList("downstream_gene_variant", "upstream_gene_variant", "intergenic_variant", "intron_variant"));

    private Organism organism = null;
    private DataSet dataSet = null;
    private DataSource dataSource = null;
    private Ontology ontology = null;
    private String assemblyVersion = "ARS-UCD1.2";
    private String geneSource = null;

    HashSet < Transcript > transcriptSet = new HashSet < Transcript > ();
    HashSet < Consequence > consequenceSet = new HashSet < Consequence > ();
    HashSet < String > previousGeneSet = new HashSet < String > ();
    HashSet < String > seenSet = new HashSet < String > ();

    private Map < String, Consequence > consequences = new HashMap < String, Consequence > ();
    private Map < Transcript, HashSet < SequenceAlteration >> transcriptToSequenceAlterationMap = new HashMap < Transcript, HashSet < SequenceAlteration >> ();
    private Map < String, SOTerm > createdSotermMap = new HashMap < String, SOTerm > ();
    private Map < String, Gene > createdGeneMap = new HashMap < String, Gene > ();
    private Map < String, Chromosome > createdChromosomeMap = new HashMap < String, Chromosome > ();
    private Map < String, Transcript > createdTranscriptMap = new HashMap < String, Transcript > ();
    private Map < String, ConsequenceType > consequenceTypeMap = new HashMap < String, ConsequenceType > ();
    private HashMap < Integer, InterMineObject > imoTracker = new HashMap < Integer, InterMineObject > ();

    /**
     * Set the assembly version for chromosomes
     * @param assemblyVersion
     */
    public void setAssemblyVersion(String assemblyVersion) {
        this.assemblyVersion = assemblyVersion;
    }


    public void setGeneSource(String geneSource) {
        this.geneSource = geneSource;
    }

    public String getGeneSource() {
        return this.geneSource;
    }


    /**
     * Setter for taxonId (used when setting taxonId from project.xml)
     * @param taxonId
     */
    public void setTaxonId(String taxonId) {
        try {
            taxonId = String.valueOf(taxonId);
        } catch (NumberFormatException e) {
            System.out.println("Invalid taxonId: " + taxonId);
        }
    }

    /**
     * Called by parent process() method for each file found
     *
     * {@inheritDoc}
     */
    public void processFile(File theFile) {
        System.out.println("Processing file: " + theFile.getName());
        if (!theFile.getName().endsWith(".vcf")) {
            LOG.info("Ignoring file " + theFile.getName() + ". Not a VCF file");
        } else {
            if (taxonId == null) {
                throw new BuildException("taxonId must be set");
            }

            //pre-populates proxy references for all Genes, Transcripts and Chromosomes in the production database
            //preFill(genes, Gene.class);
            //preFill(transcripts, Transcript.class);
            //preFill(chromosomes, Chromosome.class);

            Iterator < String[] > lineIter = null;
            try {
                FileReader reader = new FileReader(theFile);
                // Note: FormattedTextParser will not return any of the VCF headers
                lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
            } catch (Exception e) {
                throw new BuildException("Cannot parse file: " + theFile.getName() + "\n" + e);
            }

            long counter = 0;
            while (lineIter.hasNext()) {
                String[] line = lineIter.next();
                counter++;
                try {
                    if (!processVcfEntry(line)) {
                        System.out.println("Error processing line:");
                        printLine(line);
                        return;
                    }
                } catch (ObjectStoreException e) {
                    throw new BuildException("Error processing data: " + e);
                }
                if ((counter % 100000) == 0) {
                    System.out.println("Processed " + counter + " lines...");
                    LOG.info("Processed " + counter + " lines...");
                }
            }

            // storing all delayed entities
            if (transcriptToSequenceAlterationMap.size() > 0) {
                System.out.println("Doing bulk store");
                try {
                    bulkStore();
                } catch (Exception e) {
                    System.out.println(e);
                }
            }

            if (imoTracker.size() > 0) {
                // sanity check
                System.out.println("Number of entities left behind: " + imoTracker.size());
                System.out.println("Check logs");
                for (Map.Entry < Integer, InterMineObject > entry: imoTracker.entrySet()) {
                    LOG.info("ID: " + entry.getKey() + " <=> " + entry.getValue());
                }
            }
        }
    }

    /**
     * Process VCF entry
     * @param fields
     * @return boolean
     * @throws ObjectStoreException
     */
    private boolean processVcfEntry(String[] fields) throws ObjectStoreException {
        if (fields.length < expectedHeaders.length) {
            System.out.println("Fields length of " + fields.length + " vs Expected length of " + expectedHeaders.length);
            throw new BuildException("Unexpected number of columns in VCF file");
        }

        String chromosomeIdentifier = fields[0];
        Integer position = new Integer(fields[1]);
        String id = fields[2];
        System.out.println("RSID: " + id);
        String ref = fields[3];
        String alt = fields[4];
        String info = fields[7];
        ArrayList < String > infoElements = splitInfoField(info);

        Chromosome chromosome = getChromosome(chromosomeIdentifier);
        try {
            chromosome.setAssembly(assemblyVersion);
        } catch (Error e) {}
        //String type = getKeyValuePair(infoElements, "VC=").get(1);
        String type = getKeyValuePair(infoElements, "TSA=").get(1);
        //String dbSnpBuild = getKeyValuePair(infoElements, "dbSNPBuildID=").get(1);
        String dbSnpBuild = infoElements.get(0);
        ArrayList < String > aliases = getKeyValuePair(infoElements, "alias=");
        //ArrayList<String> variantEffect = getKeyValuePair(infoElements, "VE=");
        //ArrayList<String> sift = getKeyValuePair(infoElements, "Sift=");
        ArrayList < String > csq = getKeyValuePair(infoElements, "CSQ=");
        //ArrayList<String> varPep = getKeyValuePair(infoElements, "VarPep=");
        //ArrayList<String> refPep = getKeyValuePair(infoElements, "RefPep=");
        ArrayList < HashSet > returnObject = new ArrayList < HashSet > ();
        HashSet < Consequence > consequenceSet = new HashSet < Consequence > ();
        HashSet < Transcript > transcriptSet = new HashSet < Transcript > ();

        boolean storeFeature = true;

        if (type.toUpperCase().equals("SNV")) {
            SNP snp = getDirectDataLoader().createObject(SNP.class);
            snp.setSequenceOntologyTerm(getSoTerm("SNP"));
            snp.setName("SNP");
            imoTracker.put(snp.getId(), snp);
            snp.setOrganism(getOrganism());
            snp.setReferenceAllele(ref);
            snp.setLength(1);
            snp.setAlternateAllele(alt);
            snp.setPrimaryIdentifier(id);
            snp.setChromosome(chromosome);
            // set SNP -> dataSets collection
            snp.setDataSets(new HashSet < DataSet > (Arrays.asList(getDataSet())));

            // Processing ssIds and probeIds corresponding to each rsId
            Set < AliasName > ssIdSet = null;

            if (aliases != null) {
                ssIdSet = processDbsnpSsIdentifiers((SequenceAlteration) snp, aliases.get(1));
            }

            Set < AliasName > setOfAlias = new HashSet < AliasName > ();
            if (ssIdSet != null) {
                setOfAlias.addAll(ssIdSet);
            }

            // set SNP -> alias collection
            if (setOfAlias.size() > 0) {
                snp.setAliases(setOfAlias);
            }

            // All variants in the dbSNP VCF are w.r.t. the + strand
            Location location = createLocation(chromosome, snp, position, ref, alt, "1");
            snp.setChromosomeLocation(location);

            if(csq != null){
                returnObject = processCSQ((SequenceAlteration) snp, csq.get(1));
                consequenceSet = returnObject.get(0);
                transcriptSet = returnObject.get(1);
                if (consequenceSet != null) {
                    // set SNP -> consequences collection
                    snp.setConsequences(consequenceSet);
                }

                if (transcriptSet != null) {
                    // set SNP -> transcripts collection
                    snp.setTranscripts(transcriptSet);
                    for (Transcript transcript: transcriptSet) {
                        // keeping track of all SNPs that have altered a transcript
                        if (transcriptToSequenceAlterationMap.containsKey(transcript)) {
                            // if there is even one transcript from a previous entry
                            // do not attempt a store for current SNP
                            storeFeature = false;
                            transcriptToSequenceAlterationMap.get(transcript).add(snp);
                        } else {
                            transcriptToSequenceAlterationMap.put(transcript, new HashSet < SequenceAlteration > (Arrays.asList(snp)));
                            storeFeature = false;
                        }
                    }
                }
            }

            if (storeFeature) {
                try {
                    // store the SNP
                    getDirectDataLoader().store(snp);
                    imoTracker.remove(snp.getId());
                } catch (ObjectStoreException e) {
                    throw new BuildException("Problem storing SNP: " + e);
                }
            }
        } else if (type.toUpperCase().equals("INSERTION") || type.toUpperCase().equals("DELETION") || type.toUpperCase().equals("INDEL")) {
            Indel indel = getDirectDataLoader().createObject(Indel.class);
            indel.setSequenceOntologyTerm(getSoTerm("indel"));
            indel.setName("INDEL");
            imoTracker.put(indel.getId(), indel);
            indel.setOrganism(getOrganism());
            indel.setReferenceAllele(ref);
            indel.setAlternateAllele(alt);
            indel.setPrimaryIdentifier(id);
            indel.setChromosome(chromosome);
            // set INDEL -> dataSet collection
            indel.setDataSets(new HashSet < DataSet > (Arrays.asList(getDataSet())));

            // Processing ssId and probeId corresponding to each rsId
            Set < AliasName > ssIdSet = null;

            if (aliases != null) {
                ssIdSet = processDbsnpSsIdentifiers((SequenceAlteration) indel, aliases.get(1));
            }

            Set < AliasName > setOfAlias = new HashSet < AliasName > ();
            if (ssIdSet != null) {
                setOfAlias.addAll(ssIdSet);
            }

            if (setOfAlias.size() > 0) {
                // set INDEL -> alias collection
                indel.setAliases(setOfAlias);
            }

            // All variants in the dbSNP VCF are w.r.t. the + strand
            Location location = createLocation(chromosome, indel, position, ref, alt, "1");
            indel.setChromosomeLocation(location);

            // Processing variant annotation for each rsId
            if (csq != null) {
                returnObject = processCSQ((SequenceAlteration) indel, csq.get(1));
                consequenceSet = returnObject.get(0);
                transcriptSet = returnObject.get(1);

                if (consequenceSet != null) {
                    // set INDEL -> consequences collection
                    indel.setConsequences(consequenceSet);
                }

                if (transcriptSet != null) {
                    // set Indel -> transcripts collection
                    indel.setTranscripts(transcriptSet);
                    for (Transcript transcript: transcriptSet) {
                        // keeping track of all indels that have altered a transcript
                        if (transcriptToSequenceAlterationMap.containsKey(transcript)) {
                            // if there is even one transcript from a previous entry
                            // do not attempt a store for current indel
                            storeFeature = false;
                            transcriptToSequenceAlterationMap.get(transcript).add(indel);
                        } else {
                            transcriptToSequenceAlterationMap.put(transcript, new HashSet < SequenceAlteration > (Arrays.asList(indel)));
                            storeFeature = false;
                        }
                    }
                }
            }

            if (storeFeature) {
                try {
                    // store the indel
                    getDirectDataLoader().store(indel);
                    imoTracker.remove(indel.getId());
                } catch (ObjectStoreException e) {
                    throw new BuildException("Problem storing indel: " + e);
                }
            }
        } else if (type.toUpperCase().equals("SEQUENCE_ALTERATION")) {

            SequenceAlteration tsa = getDirectDataLoader().createObject(SequenceAlteration.class);
            tsa.setSequenceOntologyTerm(getSoTerm("sequence_alteration"));
            tsa.setName("SEQUENCE_ALTERATION");
            imoTracker.put(tsa.getId(), tsa);
            tsa.setOrganism(getOrganism());
            tsa.setReferenceAllele(ref);
            tsa.setAlternateAllele(alt);
            tsa.setPrimaryIdentifier(id);
            tsa.setChromosome(chromosome);
            // set tsa -> dataSet collection
            tsa.setDataSets(new HashSet < DataSet > (Arrays.asList(getDataSet())));

            // Processing ssId and probeId corresponding to each rsId
            Set < AliasName > ssIdSet = null;

            if (aliases != null) {
                ssIdSet = processDbsnpSsIdentifiers((SequenceAlteration) tsa, aliases.get(1));
            }

            Set < AliasName > setOfAlias = new HashSet < AliasName > ();
            if (ssIdSet != null) {
                setOfAlias.addAll(ssIdSet);
            }

            if (setOfAlias.size() > 0) {
                // set tsa -> alias collection
                tsa.setAliases(setOfAlias);
            }

            // All variants in the dbSNP VCF are w.r.t. the + strand
            Location location = createLocation(chromosome, tsa, position, ref, alt, "1");
            tsa.setChromosomeLocation(location);

            // Processing variant annotation for each rsId
            if (csq != null) {
                returnObject = processCSQ((SequenceAlteration) tsa, csq.get(1));
                // returnObject = processVariantEffects((SequenceAlteration) tsa, variantEffect.get(1), csqStr, siftStr, varPepStr, refPepStr);
                consequenceSet = returnObject.get(0);
                transcriptSet = returnObject.get(1);

                if (consequenceSet != null) {
                    // set tsa -> consequences collection
                    tsa.setConsequences(consequenceSet);
                }

                if (transcriptSet != null) {
                    // set tsa -> transcripts collection
                    tsa.setTranscripts(transcriptSet);
                    for (Transcript transcript: transcriptSet) {
                        // keeping track of all tsas that have altered a transcript
                        if (transcriptToSequenceAlterationMap.containsKey(transcript)) {
                            // if there is even one transcript from a previous entry
                            // do not attempt a store for current tsa
                            storeFeature = false;
                            transcriptToSequenceAlterationMap.get(transcript).add(tsa);
                        } else {
                            transcriptToSequenceAlterationMap.put(transcript, new HashSet < SequenceAlteration > (Arrays.asList(tsa)));
                            storeFeature = false;
                        }
                    }
                }
            }

            if (storeFeature) {
                try {
                    // store the tsa
                    getDirectDataLoader().store(tsa);
                    imoTracker.remove(tsa.getId());
                } catch (ObjectStoreException e) {
                    throw new BuildException("Problem storing Tandem Repeat: " + e);
                }
            }
        } else if (type.toUpperCase().equals("SUBSTITUTION")) {
            Substitution tsa = getDirectDataLoader().createObject(Substitution.class);
            tsa.setSequenceOntologyTerm(getSoTerm("substitution"));
            tsa.setName("SUBSTITUTION");
            imoTracker.put(tsa.getId(), tsa);
            tsa.setOrganism(getOrganism());
            tsa.setReferenceAllele(ref);
            tsa.setAlternateAllele(alt);
            tsa.setPrimaryIdentifier(id);
            tsa.setChromosome(chromosome);
            // set tsa -> dataSet collection
            tsa.setDataSets(new HashSet < DataSet > (Arrays.asList(getDataSet())));

            // Processing ssId and probeId corresponding to each rsId
            Set < AliasName > ssIdSet = null;

            if (aliases != null) {
                ssIdSet = processDbsnpSsIdentifiers((SequenceAlteration) tsa, aliases.get(1));
            }

            Set < AliasName > setOfAlias = new HashSet < AliasName > ();
            if (ssIdSet != null) {
                setOfAlias.addAll(ssIdSet);
            }

            if (setOfAlias.size() > 0) {
                // set tsa -> alias collection
                tsa.setAliases(setOfAlias);
            }

            // All variants in the dbSNP VCF are w.r.t. the + strand
            Location location = createLocation(chromosome, tsa, position, ref, alt, "1");
            tsa.setChromosomeLocation(location);

            // Processing variant annotation for each rsId
            if (csq != null) {
                returnObject = processCSQ((SequenceAlteration) tsa, csq.get(1));
                consequenceSet = returnObject.get(0);
                transcriptSet = returnObject.get(1);

                if (consequenceSet != null) {
                    // set tsa -> consequences collection
                    tsa.setConsequences(consequenceSet);
                }

                if (transcriptSet != null) {
                    // set tsa -> transcripts collection
                    tsa.setTranscripts(transcriptSet);
                    for (Transcript transcript: transcriptSet) {
                        // keeping track of all tsas that have altered a transcript
                        if (transcriptToSequenceAlterationMap.containsKey(transcript)) {
                            // if there is even one transcript from a previous entry
                            // do not attempt a store for current tsa
                            storeFeature = false;
                            transcriptToSequenceAlterationMap.get(transcript).add(tsa);
                        } else {
                            transcriptToSequenceAlterationMap.put(transcript, new HashSet < SequenceAlteration > (Arrays.asList(tsa)));
                            storeFeature = false;
                        }
                    }
                }
            }

            if (storeFeature) {
                try {
                    // store the tsa
                    getDirectDataLoader().store(tsa);
                    imoTracker.remove(tsa.getId());
                } catch (ObjectStoreException e) {
                    throw new BuildException("Problem storing Tandem Repeat: " + e);
                }
            }
        } else if (type.toUpperCase().equals("TANDEM_REPEAT")) {


            TandemRepeat tsa = getDirectDataLoader().createObject(TandemRepeat.class);
            tsa.setSequenceOntologyTerm(getSoTerm("tandem_repeat"));
            tsa.setName("TANDEM_REPEAT");
            imoTracker.put(tsa.getId(), tsa);
            tsa.setOrganism(getOrganism());
            tsa.setReferenceAllele(ref);
            tsa.setAlternateAllele(alt);
            tsa.setPrimaryIdentifier(id);
            tsa.setChromosome(chromosome);
            // set tsa -> dataSet collection
            tsa.setDataSets(new HashSet < DataSet > (Arrays.asList(getDataSet())));

            // Processing ssId and probeId corresponding to each rsId
            Set < AliasName > ssIdSet = null;

            if (aliases != null) {
                ssIdSet = processDbsnpSsIdentifiers((SequenceAlteration) tsa, aliases.get(1));
            }

            Set < AliasName > setOfAlias = new HashSet < AliasName > ();
            if (ssIdSet != null) {
                setOfAlias.addAll(ssIdSet);
            }

            if (setOfAlias.size() > 0) {
                // set tsa -> alias collection
                tsa.setAliases(setOfAlias);
            }

            // All variants in the dbSNP VCF are w.r.t. the + strand
            Location location = createLocation(chromosome, tsa, position, ref, alt, "1");
            tsa.setChromosomeLocation(location);

            // Processing variant annotation for each rsId
            if (csq != null) {
                returnObject = processCSQ((SequenceAlteration) tsa, csq.get(1));
                consequenceSet = returnObject.get(0);
                transcriptSet = returnObject.get(1);

                if (consequenceSet != null) {
                    // set tsa -> consequences collection
                    tsa.setConsequences(consequenceSet);
                }

                if (transcriptSet != null) {
                    // set tsa -> transcripts collection
                    tsa.setTranscripts(transcriptSet);
                    for (Transcript transcript: transcriptSet) {
                        // keeping track of all tsas that have altered a transcript
                        if (transcriptToSequenceAlterationMap.containsKey(transcript)) {
                            // if there is even one transcript from a previous entry
                            // do not attempt a store for current tsa
                            storeFeature = false;
                            transcriptToSequenceAlterationMap.get(transcript).add(tsa);
                        } else {
                            transcriptToSequenceAlterationMap.put(transcript, new HashSet < SequenceAlteration > (Arrays.asList(tsa)));
                            storeFeature = false;
                        }
                    }
                }
            }

            if (storeFeature) {
                try {
                    // store the tsa
                    getDirectDataLoader().store(tsa);
                    imoTracker.remove(tsa.getId());
                } catch (ObjectStoreException e) {
                    throw new BuildException("Problem storing Tandem Repeat: " + e);
                }
            }
        } else {
            System.out.println("Unexpected feature of type: " + type);
            System.exit(1);
        }
        return true;
    }

    /**
     * Stores all delayed entities
     * @throws ObjectStoreException
     */
    private void bulkStore() throws ObjectStoreException {
        for (Map.Entry < Transcript, HashSet < SequenceAlteration >> entry: transcriptToSequenceAlterationMap.entrySet()) {
            Transcript transcript = entry.getKey();
            // set Transcript -> sequenceAlterations collection
            transcript.setSequenceAlterations(entry.getValue());

            for (SequenceAlteration sa: entry.getValue()) {
                if (seenSet.contains(String.valueOf(sa.getId()))) {
                    if (imoTracker.containsKey(sa.getId())) {
                        System.out.println("imoTracker has sequenceAlteration even after its been stored: " + sa);
                    }
                    continue;
                } else {
                    // store sequence alteration since it was never seen before
                    getDirectDataLoader().store(sa);
                    imoTracker.remove(sa.getId());
                    seenSet.add(String.valueOf(sa.getId()));
                }
            }

            getDirectDataLoader().store(transcript);
            imoTracker.remove(transcript.getId());
        }
        transcriptToSequenceAlterationMap.clear();
    }

    /**
     * Creates a set of AliasName entities for a given series of ssIds
     * @param saFeature
     * @param ssIds
     * @return set of <AliasName>
     * @throws ObjectStoreException
     */
    private Set < AliasName > processDbsnpSsIdentifiers(SequenceAlteration saFeature, String ssIds) throws ObjectStoreException {
        // Note: The genomic_model will have to be altered where reference to AliasName is now in BioEntity instead of Gene
        Set < AliasName > setOfSsIdObjects = new HashSet < AliasName > ();
        String[] aliases = ssIds.split(",");
        for (int i = 0; i < aliases.length; i++) {
            String[] eachAlias = aliases[i].split(":");
            String ssId = eachAlias[0];
            String ssIdSource = eachAlias[1];
            AliasName alias = getDirectDataLoader().createObject(AliasName.class);
            imoTracker.put(alias.getId(), alias);
            alias.setIdentifier(ssId);
            alias.setSource(ssIdSource);
            alias.setOrganism(getOrganism());
            // set AliasName -> feature reference
            alias.setFeatures(new HashSet < SequenceFeature > (Arrays.asList((SequenceFeature) saFeature)));
            getDirectDataLoader().store(alias);
            imoTracker.remove(alias.getId());
            setOfSsIdObjects.add(alias);
        }
        return setOfSsIdObjects;
    }

    /**
     * Creates a set of AliasName entities for a given series of probeIds
     * @param saFeature
     * @param probeIds
     * @return set of <AliasName>
     * @throws ObjectStoreException
     */
    private Set < AliasName > processSnpArrayProbeIdentifiers(SequenceAlteration saFeature, String probeIds) throws ObjectStoreException {
        Set < AliasName > setOfProbeIdObjects = new HashSet < AliasName > ();
        String[] probes = probeIds.split(",");
        for (int i = 0; i < probes.length; i++) {
            String[] eachProbe = probes[i].split(":");
            String probeId = eachProbe[0];
            String probeSource = eachProbe[1];
            AliasName alias = getDirectDataLoader().createObject(AliasName.class);
            imoTracker.put(alias.getId(), alias);
            alias.setIdentifier(probeId);
            alias.setSource(probeSource);
            alias.setOrganism(getOrganism());
            // AliasName -> feature reference
            alias.setFeatures(new HashSet < SequenceFeature > (Arrays.asList((SequenceFeature) saFeature)));
            getDirectDataLoader().store(alias);
            imoTracker.remove(alias.getId());
            setOfProbeIdObjects.add(alias);
        }
        return setOfProbeIdObjects;
    }

    /**
     * Processes variant effects
     * @param saFeature
     * @param variantEffects
     * @return
     * @throws ObjectStoreException
     * returnObject = processVariantEffects((SequenceAlteration) snp, variantEffect.get(1), csqStr, siftStr, varPepStr, refPepStr);
     */
    private ArrayList < HashSet > processCSQ(SequenceAlteration saFeature, String csq) throws ObjectStoreException {
        String[] csqs = csq.split(",");

        HashSet < Consequence > consequenceSet = new HashSet < Consequence > ();
        HashSet < Transcript > transcriptSet = new HashSet < Transcript > ();

        for (int i = 0; i < csqs.length; i++) {
            // order of info: [ GeneID, TranscriptID, consequence, impact, cDNA_pos, CDS_pos, Protein_pos, AA_change, Codon_change ]
            HashSet < ConsequenceType > consequenceTypeSet = new HashSet < ConsequenceType > ();
            String[] csqInfo = csqs[i].split("\\|");
		System.out.println("CSQ: " + Arrays.toString(csqInfo));
            /*
                     if (FUNCTION_CLASS_TO_IGNORE.contains(annotationInfo[2])) {
                         continue;
                     }
                     else {
                         String[] consequenceTypes = annotationInfo[2].split("\\|");
                         for (String consequenceType : consequenceTypes) {
                             consequenceTypeSet.add(getConsequence(consequenceType));
                         }
                     }
                 */
            //consequence.setCdnaPosition(annotationInfo[4]);
            //consequence.setCdsPosition(annotationInfo[5]);
            //consequence.setAminoAcidPosition(annotationInfo[6]);

            // set Consequence -> consequenceType collection
            // TODO: Issue on the webapp while trying to access consequenceType
            Consequence consequence = getDirectDataLoader().createObject(Consequence.class);
            imoTracker.put(consequence.getId(), consequence);
            consequence.setRsId(saFeature.getPrimaryIdentifier());
            consequence.setConsequenceTypes(consequenceTypeSet);

            if (csqInfo.length > 1) {

                consequence.setAlternateCodon(csqInfo[0]);
                consequence.setAlternateAllele(csqInfo[0]);
                consequenceTypeSet.add(getConsequence(csqInfo[1]));

                //System.out.println("csqs INFOR: " + csqInfo.length + " " + csqInfo.toString());
                if(csqInfo.length > 3) {
                	consequence.setTranscriptIdentifier(csqInfo[3]);
               	 	Transcript transcript = getTranscript(csqInfo[3]);
                	consequence.setTranscript(transcript);
                	transcriptSet.add(transcript);
		}

                if(csqInfo.length > 4 && csqInfo[4] != null){
                    String[] residue = csqInfo[4].split("\\/");
                    consequence.setReferenceResidue(residue[0]);
                    consequence.setAlternateResidue(residue[1]);
                }

                if(csqInfo.length > 5){
                    String[] sift = csqInfo[5].split("\\(");
                    //sift[1] = sift[1].split("\\)", "");
                    consequence.setSiftQualitativePrediction(sift[0].replaceAll("-_", ""));
                    consequence.setSiftNumericalValue(sift[1].replaceAll("[()]", ""));
                }

                

                // if (sift != null) {
                //     if (siftArr.length > 0 && i < siftArr.length) {
                //         String[] siftSplit = siftArr[i].split("\\|");
                //         consequence.setSiftQualitativePrediction(siftSplit[1]);
                //         consequence.setSiftNumericalValue(siftSplit[2]);
                //     } else {
                //         String[] siftSplit = sift.split("\\|");
                //         consequence.setSiftQualitativePrediction(siftSplit[1]);
                //         consequence.setSiftNumericalValue(siftSplit[2]);
                //     }
                // }

                // if (varPep != null) {
                //     if (varPepArr.length > 0 && i < varPepArr.length) {
                //         String[] varPepSplit = varPepArr[i].split("\\|");
                //         consequence.setAlternateResidue(varPepSplit[1]);
                //     } else {
                //         String[] varPepSplit = varPep.split("\\|");
                //         consequence.setAlternateResidue(varPepSplit[1]);
                //     }
                // }

            } 

            // if (refPep != null) {
            //     consequence.setReferenceResidue(refPep);
            // }
            /*
                     String[] substitution = annotationInfo[7].split("/");
                     if (substitution.length > 0) {
                         consequence.setReferenceResidue(substitution[0]);
                         if (substitution.length != 2) {
                             LOG.warn("Consequence type " + annotationInfo[2] + " with just one residue info: " + annotationInfo[7] + ". Treating it as referenceResidue");
                         } else {
                             consequence.setAlternateResidue(substitution[1]);
                         }
                     }

                     if (! "-".equals(annotationInfo[8])) {
                         String[] codons = annotationInfo[8].split("/");
                         if (codons.length > 0) {
                             consequence.setReferenceCodon(codons[0].toUpperCase());
                             consequence.setAlternateCodon(codons[1].toUpperCase());
                         }
                     }
                 */
            // set Consequence -> variant reference
            consequence.setVariant(saFeature);
            getDirectDataLoader().store(consequence);
            imoTracker.remove(consequence.getId());
            consequenceSet.add(consequence);
        }
        return new ArrayList < HashSet > (Arrays.asList(consequenceSet, transcriptSet));
    }

    /**
     * Processes variant effects
     * @param saFeature
     * @param variantEffects
     * @return
     * @throws ObjectStoreException
     * returnObject = processVariantEffects((SequenceAlteration) snp, variantEffect.get(1), csqStr, siftStr, varPepStr, refPepStr);
     */
    private ArrayList < HashSet > processVariantEffects(SequenceAlteration saFeature, String variantEffects, String csq, String sift, String varPep, String refPep) throws ObjectStoreException {
        String[] effects = variantEffects.split(",");
        String[] siftArr = null;
        String[] csqArr = null;
        String[] varPepArr = null;

        if (sift != null) {
            siftArr = sift.split(",");
        }
        if (csq != null) {
            csqArr = csq.split(",");
        }
        if (varPep != null) {
            varPepArr = varPep.split(",");
        }

        HashSet < Consequence > consequenceSet = new HashSet < Consequence > ();
        HashSet < Transcript > transcriptSet = new HashSet < Transcript > ();

        for (int i = 0; i < effects.length; i++) {
            // order of info: [ GeneID, TranscriptID, consequence, impact, cDNA_pos, CDS_pos, Protein_pos, AA_change, Codon_change ]
            HashSet < ConsequenceType > consequenceTypeSet = new HashSet < ConsequenceType > ();
            String[] effectInfo = effects[i].split("\\|");
            /*
                     if (FUNCTION_CLASS_TO_IGNORE.contains(annotationInfo[2])) {
                         continue;
                     }
                     else {
                         String[] consequenceTypes = annotationInfo[2].split("\\|");
                         for (String consequenceType : consequenceTypes) {
                             consequenceTypeSet.add(getConsequence(consequenceType));
                         }
                     }
                 */
            //consequence.setCdnaPosition(annotationInfo[4]);
            //consequence.setCdsPosition(annotationInfo[5]);
            //consequence.setAminoAcidPosition(annotationInfo[6]);

            // set Consequence -> consequenceType collection
            // TODO: Issue on the webapp while trying to access consequenceType
            Consequence consequence = getDirectDataLoader().createObject(Consequence.class);
            imoTracker.put(consequence.getId(), consequence);
            consequence.setRsId(saFeature.getPrimaryIdentifier());
            consequence.setConsequenceTypes(consequenceTypeSet);

            if (effectInfo.length > 1) {
                //System.out.println("EFFECTS INFOR: " + effectInfo.length + " " + effectInfo.toString());
                consequenceTypeSet.add(getConsequence(effectInfo[0]));
                Transcript transcript = getTranscript(effectInfo[3]);
                consequence.setTranscriptIdentifier(effectInfo[3]);
                consequence.setTranscript(transcript);
                transcriptSet.add(transcript);

                if (sift != null) {
                    if (siftArr.length > 0 && i < siftArr.length) {
                        String[] siftSplit = siftArr[i].split("\\|");
                        consequence.setSiftQualitativePrediction(siftSplit[1]);
                        consequence.setSiftNumericalValue(siftSplit[2]);
                    } else {
                        String[] siftSplit = sift.split("\\|");
                        consequence.setSiftQualitativePrediction(siftSplit[1]);
                        consequence.setSiftNumericalValue(siftSplit[2]);
                    }
                }

                if (varPep != null) {
                    if (varPepArr.length > 0 && i < varPepArr.length) {
                        String[] varPepSplit = varPepArr[i].split("\\|");
                        consequence.setAlternateResidue(varPepSplit[1]);
                    } else {
                        String[] varPepSplit = varPep.split("\\|");
                        consequence.setAlternateResidue(varPepSplit[1]);
                    }
                }

            } else {
                consequenceTypeSet.add(getConsequence(effects[i]));
            }

            if (refPep != null) {
                consequence.setReferenceResidue(refPep);
            }
            /*
                     String[] substitution = annotationInfo[7].split("/");
                     if (substitution.length > 0) {
                         consequence.setReferenceResidue(substitution[0]);
                         if (substitution.length != 2) {
                             LOG.warn("Consequence type " + annotationInfo[2] + " with just one residue info: " + annotationInfo[7] + ". Treating it as referenceResidue");
                         } else {
                             consequence.setAlternateResidue(substitution[1]);
                         }
                     }

                     if (! "-".equals(annotationInfo[8])) {
                         String[] codons = annotationInfo[8].split("/");
                         if (codons.length > 0) {
                             consequence.setReferenceCodon(codons[0].toUpperCase());
                             consequence.setAlternateCodon(codons[1].toUpperCase());
                         }
                     }
                 */
            // set Consequence -> variant reference
            consequence.setVariant(saFeature);
            getDirectDataLoader().store(consequence);
            imoTracker.remove(consequence.getId());
            consequenceSet.add(consequence);
        }
        return new ArrayList < HashSet > (Arrays.asList(consequenceSet, transcriptSet));
    }

    /**
     * Creates a Location entity for a given set of location information
     * @param locatedOn
     * @param feature
     * @param start
     * @param refAllele
     * @param altAllele
     * @param strand
     * @throws ObjectStoreException
     */
    private Location createLocation(Chromosome locatedOn, SequenceAlteration feature, int start, String refAllele, String altAllele, String strand) throws ObjectStoreException {
        int length = 0;
        Location location = getDirectDataLoader().createObject(Location.class);
        imoTracker.put(location.getId(), location);
        location.setLocatedOn(locatedOn);
        location.setFeature((BioEntity) feature);
        length = refAllele.length();

        int end = (start + length) - 1;

        if (start <= end) {
            location.setStart(start);
            location.setEnd(end);
        } else {
            System.out.println("Trying to create a Location entity where start > end");
            System.exit(1);
        }

        location.setStrand(strand);
        // Annotating Location entity such that it is not considered when computing overlaps during post-process
        location.setDoNotComputeOverlaps("Y");
        getDirectDataLoader().store(location);
        imoTracker.remove(location.getId());
        return location;
    }

    /**
     * Splits INFO fields of a VCF entry
     * @param infoField
     * @return
     */
    private ArrayList < String > splitInfoField(String infoField) {
        return new ArrayList < String > (Arrays.asList(StringUtil.split(infoField, ";")));
    }

    /**
     * For a given tag, returns its corresponding value
     * @param elements
     * @param tag
     * @return
     */
    private ArrayList < String > getKeyValuePair(ArrayList < String > elements, String tag) {
        for (int i = 0; i < elements.size(); i++) {
            if (elements.get(i).contains(tag)) {
                return new ArrayList < String > (Arrays.asList(StringUtil.split(elements.get(i), "=")));
            }
        }
        return null;
    }

    /**
     * Returns the current DataSource
     * @return
     * @throws ObjectStoreException
     */
    private DataSource getDataSource() throws ObjectStoreException {
        if (dataSource == null) {
            dataSource = getDirectDataLoader().createObject(DataSource.class);
            imoTracker.put(dataSource.getId(), dataSource);
            dataSource.setName(DATA_SOURCE_NAME);
            getDirectDataLoader().store(dataSource);
            imoTracker.remove(dataSource.getId());
        }
        return dataSource;
    }

    /**
     * Returns the current DataSet
     * @return
     * @throws ObjectStoreException
     */
    private DataSet getDataSet() throws ObjectStoreException {
        if (dataSet == null) {
            dataSet = getDirectDataLoader().createObject(DataSet.class);
            imoTracker.put(dataSet.getId(), dataSet);
            dataSet.setName(DATASET_TITLE);
            dataSet.setDataSource(getDataSource());
            getDirectDataLoader().store(dataSet);
            imoTracker.remove(dataSet.getId());
        }
        return dataSet;
    }

    /**
     * Returns the current Organism
     * @return organism
     * @throws ObjectStoreException
     */
    private Organism getOrganism() throws ObjectStoreException {
        if (organism == null) {
            organism = getDirectDataLoader().createObject(Organism.class);
            imoTracker.put(organism.getId(), organism);
            organism.setTaxonId(taxonId);
            getDirectDataLoader().store(organism);
            imoTracker.remove(organism.getId());
        }
        return organism;
    }

    /**
     * Returns a ConsequenceType entity for a given consequence type
     * @param consequenceTypeString
     * @return
     * @throws ObjectStoreException
     */
    private ConsequenceType getConsequence(String consequenceTypeString) throws ObjectStoreException {
        ConsequenceType consequenceType;
        if (consequenceTypeMap.containsKey(consequenceTypeString)) {
            consequenceType = consequenceTypeMap.get(consequenceTypeString);
        } else {
            consequenceType = getDirectDataLoader().createObject(ConsequenceType.class);
            consequenceType.setName(consequenceTypeString.toUpperCase());
            getDirectDataLoader().store(consequenceType);
            consequenceTypeMap.put(consequenceTypeString, consequenceType);
        }
        return consequenceType;
    }

    /**
     * For a given identifier, returns a Gene entity
     * @param identifier
     * @return
     * @throws ObjectStoreException
     */
    private Gene getGene(String identifier) throws ObjectStoreException {
        Gene gene;
        if (createdGeneMap.containsKey(identifier)) {
            gene = createdGeneMap.get(identifier);
        } else {
            gene = getDirectDataLoader().createObject(Gene.class);
            gene.setSequenceOntologyTerm(getSoTerm("gene"));
            imoTracker.put(gene.getId(), gene);
            // TODO: Might have to set source of Gene as well, to enable integration downstream
            gene.setOrganism(getOrganism());
            createdGeneMap.put(identifier, gene);
            //getDirectDataLoader().store(gene);
        }
        return gene;
    }

    /**
     * For a given identifier, returns a Transcript entity
     * @param identifier
     * @return
     * @throws ObjectStoreException
     */
    private Transcript getTranscript(String identifier) throws ObjectStoreException {
        Transcript transcript;
        if (createdTranscriptMap.containsKey(identifier)) {
            transcript = createdTranscriptMap.get(identifier);
        } else {
            transcript = getDirectDataLoader().createObject(Transcript.class);
            transcript.setSequenceOntologyTerm(getSoTerm("transcript"));
            transcript.setPrimaryIdentifier(identifier);
            transcript.setSource(this.geneSource);
            transcript.setOrganism(getOrganism());
            imoTracker.put(transcript.getId(), transcript);
            createdTranscriptMap.put(identifier, transcript);
        }
        return transcript;
    }

    /**
     * For a given identifier, returns a Chromosome entity
     * @param identifier
     * @return
     * @throws ObjectStoreException
     */
    private Chromosome getChromosome(String identifier) throws ObjectStoreException {
        Chromosome chr;
        if (createdChromosomeMap.containsKey(identifier)) {
            chr = createdChromosomeMap.get(identifier);
        } else {
            chr = getDirectDataLoader().createObject(Chromosome.class);
            imoTracker.put(chr.getId(), chr);
            chr.setPrimaryIdentifier(identifier);
            chr.setOrganism(getOrganism());
            chr.addDataSets(getDataSet());
            chr.setSequenceOntologyTerm(getSoTerm("chromosome"));
            getDirectDataLoader().store(chr);
            imoTracker.remove(chr.getId());
            createdChromosomeMap.put(identifier, chr);
        }
        return chr;
    }

    /**
     * Returns a Ontology entity for 'Sequence Ontology'
     * @return
     * @throws ObjectStoreException
     */
    private Ontology getSequenceOntology() throws ObjectStoreException {
        if (ontology == null) {
            ontology = getDirectDataLoader().createObject(Ontology.class);
            imoTracker.put(ontology.getId(), ontology);
            ontology.setName("Sequence Ontology");
            ontology.setUrl("http://www.sequenceontology.org");
            getDirectDataLoader().store(ontology);
            imoTracker.remove(ontology.getId());
        }
        return ontology;
    }

    /**
     * Returns a SOTerm entity for a given SO feature type
     * @param featureType
     * @return
     * @throws ObjectStoreException
     */
    private SOTerm getSoTerm(String featureType) throws ObjectStoreException {
        SOTerm soTerm = createdSotermMap.get(featureType);
        if (soTerm == null) {
            System.out.print("FEATURE TYPE " + featureType);
            soTerm = getDirectDataLoader().createObject(SOTerm.class);
            imoTracker.put(soTerm.getId(), soTerm);
            soTerm.setOntology(getSequenceOntology());
            soTerm.setName(featureType);
            getDirectDataLoader().store(soTerm);
            imoTracker.remove(soTerm.getId());
            createdSotermMap.put(featureType, soTerm);
        }
        return soTerm;
    }

    /**
     * Convenience method for printing a list of string
     * @param line
     */
    private void printLine(String[] line) {
        for (int i = 0; i < line.length; i++) {
            System.out.print(line[i] + "\t");
        }
        System.out.print("\n");
    }


    /**
     * Queries the production database and returns a Map of proxy references for a entities of a given class
     * @param map
     * @param objectClass
     */
    private void preFill(Map < String, ProxyReference > map, Class < ? extends InterMineObject > objectClass) {
        Query q = new Query();
        QueryClass qC = new QueryClass(objectClass);
        q.addFrom(qC);
        QueryField qFName = new QueryField(qC, "primaryIdentifier");
        QueryField qFId = new QueryField(qC, "id");
        q.addToSelect(qFName);
        q.addToSelect(qFId);
        QueryClass qcOrg = new QueryClass(Organism.class);
        q.addFrom(qcOrg);
        QueryObjectReference orgRef = new QueryObjectReference(qC, "organism");
        QueryField qFTaxonId = new QueryField(qcOrg, "taxonId");

        ConstraintSet cs = new ConstraintSet(ConstraintOp.AND);
        cs.addConstraint(new ContainsConstraint(orgRef, ConstraintOp.CONTAINS, qcOrg));
        cs.addConstraint(new SimpleConstraint(qFTaxonId, ConstraintOp.EQUALS, new QueryValue(taxonId)));
        q.setConstraint(cs);

        LOG.info("Prefilling ProxyReferences with query: " + q);

        try {
            Results res = getIntegrationWriter().getObjectStore().execute(q, 5000, false, false, false);
            Iterator < Object > resIter = res.iterator();
            while (resIter.hasNext()) {
                ResultsRow < Object > rr = (ResultsRow < Object > ) resIter.next();
                String name = (String) rr.get(0);
                Integer id = (Integer) rr.get(1);
                map.put(name, new ProxyReference(getIntegrationWriter().getObjectStore(), id, objectClass));
            }
        } catch (Exception e) {
            throw new BuildException("Problem in prefilling ProxyReferences: " + e);
        }
        LOG.info("Prefill complete with map containing " + map.size() + " ProxyReferences");
        System.out.println("Prefill complete with map containing " + map.size() + " ProxyReferences");
    }
}
