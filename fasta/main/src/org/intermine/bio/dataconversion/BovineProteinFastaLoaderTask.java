package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2013 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */
import java.lang.IllegalAccessException;
import java.lang.Override;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.apache.tools.ant.BuildException;
import org.apache.commons.lang.StringUtils;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;
import java.util.NoSuchElementException;

import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.Sequence;
import org.intermine.metadata.Model;
import org.intermine.model.FastPathObject;
import org.intermine.model.InterMineObject;
import org.intermine.model.bio.BioEntity;
import org.intermine.model.bio.DataSet;
import org.intermine.model.bio.Location;
import org.intermine.model.bio.Organism;
import org.intermine.model.bio.SequenceFeature;
import org.intermine.objectstore.ObjectStore;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.DynamicUtil;
import org.intermine.objectstore.query.PendingClob;
import org.intermine.metadata.Util;

/**
 * A fasta loader that understand the headers of CDS fasta files and can make the
 * appropriate extra objects and references.
 * @author Kim Rutherford
 * This script has been adapted with the help of AIPCDSFastaLoaderTask.java and FlyBaseCDSFastaLoaderTask.java
 */
public class BovineProteinFastaLoaderTask extends BovineFeatureFastaLoaderTask {

    private String classAttribute = "primaryIdentifier";
    private int storeCount = 0;

    // hashmap to keep track of InterMineObject of type Gene
    private Map<String, InterMineObject> geneIdMap = new HashMap<String, InterMineObject>();
    // hashmap to keep track of InterMineObject of type MRNA
    private Map<String, InterMineObject> mrnaIdMap = new HashMap<String, InterMineObject>();

    /**
     * {@inheritDoc}
     */
    @Override
    protected void processSequence(Organism organism, Sequence bioJavaSequence) throws ObjectStoreException {

        if (organism == null) {
            return;
        }

        org.intermine.model.bio.Sequence proteinSequenceObject = getDirectDataLoader().createObject(
                org.intermine.model.bio.Sequence.class);

        String sequence = bioJavaSequence.seqString();
        String md5checksum = Util.getMd5checksum(sequence);
        proteinSequenceObject.setResidues(new PendingClob(sequence));
        proteinSequenceObject.setLength(bioJavaSequence.length());
        proteinSequenceObject.setMd5checksum(md5checksum);

        Class<? extends InterMineObject> imClass;
        Class<?> c;
        try {
            c = Class.forName(getClassName());
            if (InterMineObject.class.isAssignableFrom(c)) {
                imClass = (Class<? extends InterMineObject>) c;
            } else {
                throw new RuntimeException("Feature className must be a valid class in the model that inherits from InterMineObject, but was: " + getClassName());
            }
        } catch (ClassNotFoundException e1) {
            throw new RuntimeException("unknown class: " + getClassName() + " while creating new Sequence object");
        }

        BioEntity imo = (BioEntity) getDirectDataLoader().createObject(imClass);
        String sequenceName = getIdentifier(bioJavaSequence);
        String attributeValue = sequenceName;

        // set primaryIdentifier
        try {
            System.out.println("Setting " + classAttribute + " as " + attributeValue);
            imo.setFieldValue(classAttribute, attributeValue);
        } catch (Exception e) {
            throw new IllegalArgumentException("Error setting: " + getClassName() + "." + classAttribute + " to: " + attributeValue + ". Does the attribute exist?");
        }

        // set sequence
        try {
            imo.setFieldValue("sequence", proteinSequenceObject);
        } catch (Exception e) {
            throw new IllegalArgumentException("Error setting: " + getClassName() + ".sequence to: " + attributeValue + ". Does the attribute exist?");
        }

        // set organism
        imo.setOrganism(organism);

        // set length
        try {
            imo.setFieldValue("length", new Integer(proteinSequenceObject.getLength()));
        } catch (Exception e) {
            throw new IllegalArgumentException("Error setting: " + getClassName() + ".length to: " + proteinSequenceObject.getLength() + ". Does the attribute exist?");
        }

        // set md5checksum
        try {
            imo.setFieldValue("md5checksum", md5checksum);
        } catch (Exception e) {
        }

        extraProcessing(bioJavaSequence, proteinSequenceObject, imo, organism, getDataSet());

        DataSet dataSet = getDataSet();
        imo.addDataSets(dataSet);

        try {
            getDirectDataLoader().store(proteinSequenceObject);
            getDirectDataLoader().store(imo);
            storeCount += 2;
        } catch (ObjectStoreException e) {
            throw new BuildException("store failed", e);
        }
    }


    /**
     * {@inheritDoc}
     */
    @Override
    protected void extraProcessing(Sequence bioJavaSequence, org.intermine.model.bio.Sequence flymineSequence, BioEntity bioEntity, Organism organism, DataSet dataSet) throws ObjectStoreException {

        String geneIdentifier = null;
        String mrnaIdentifier = null;
        String source = getGeneSource();
        Annotation annotation = bioJavaSequence.getAnnotation();
        String header = (String) annotation.getProperty("description");
        String regexp = ".(\\S+)\\s(\\S+)";
        Pattern p = Pattern.compile(regexp);
        Matcher m = p.matcher(header);
        if (m.matches()) {
            geneIdentifier = m.group(1);
            mrnaIdentifier = m.group(2);
        }

        ObjectStore os = getIntegrationWriter().getObjectStore();
        Model model = os.getModel();
        // check to see if 'Polypeptide' class exists in the model
        if (model.hasClassDescriptor(model.getPackageName() + ".Polypeptide")) {
            Class<? extends FastPathObject> cdsCls = model.getClassDescriptorByName("Polypeptide").getType();
            if (!DynamicUtil.isInstance(bioEntity, cdsCls)) {
                throw new RuntimeException("the InterMineObject passed to "
                        + "RefSeqProteinFastaLoaderTask.extraProcessing() is not a "
                        + "Polypeptide: " + bioEntity);
            }
	    bioEntity.setFieldValue("source", source);
            if (geneIdentifier != null) {
                InterMineObject gene = getGene(geneIdentifier, getGeneSource(), organism, model);
                // setting the 'geneIdentifier' attribute for class 'Polypeptide'
                bioEntity.setFieldValue("geneIdentifier", geneIdentifier);
                // setting the 'gene' reference for class 'Polypeptide'
                // the reason why gene object is being passed is because setFieldValue() expects
                // an object as an argument for setting a reference
                bioEntity.setFieldValue("gene", gene);
                try {
                    // adding Polypeptide object to collection 'polypeptides' in class 'Gene'
                    HashSet polypeptidesCollection = (HashSet) gene.getFieldValue("polypeptides");
                    polypeptidesCollection.add(bioEntity);
                    gene.setFieldValue("polypeptides", polypeptidesCollection);
                    // updating the geneIdMap
                    geneIdMap.put(geneIdentifier, gene);
                } catch(IllegalAccessException e) {
                    e.printStackTrace();
                }

            }
            if (mrnaIdentifier != null) {
                InterMineObject mrna = getMRNA(mrnaIdentifier, getGeneSource(), organism, model);
                // setting the 'mrnaIdentifier' attribute for class 'Polypeptide'
                bioEntity.setFieldValue("mrnaIdentifier", mrnaIdentifier);
                // setting the 'mrna' reference for class 'Polypeptide'
                // the reason why mrna object is being passed is because setFieldValue() expects
                // an object as an argument for setting a reference
                bioEntity.setFieldValue("mrna", mrna);
                try {
                    // adding Polypeptide object to collection 'polypeptide' in class 'MRNA'
                    HashSet polypeptideCollection = (HashSet) mrna.getFieldValue("polypeptide");
                    polypeptideCollection.add(bioEntity);
                    mrna.setFieldValue("polypeptide", polypeptideCollection);
                    // updating the mrnaIdMap
                    mrnaIdMap.put(mrnaIdentifier, mrna);
                } catch(IllegalAccessException e) {
                    e.printStackTrace();
                }
            }
        } else {
            throw new RuntimeException("Trying to load Protein sequence but Protein does not exist in the data model");
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected InterMineObject getGene(String identifier, String source, Organism organism, Model model) throws ObjectStoreException {
        // overriding getGene method to get more control over how and when the Gene objects are stored
        InterMineObject gene = null;
        if (geneIdMap.containsKey(identifier)) {
            // if geneIdMap contains the given geneIdentifier as a key then fetch the gene object
            // corresponding to the key
            gene = geneIdMap.get(identifier);
        }
        else {
            // else create a new gene object and store in geneIdMap
            if (model.hasClassDescriptor(model.getPackageName() + ".Gene")) {
                @SuppressWarnings("unchecked") Class<? extends InterMineObject> geneCls = (Class<? extends InterMineObject>) model.getClassDescriptorByName("Gene").getType();
                gene = getDirectDataLoader().createObject(geneCls);
                gene.setFieldValue("primaryIdentifier", identifier);
                gene.setFieldValue("source", source);
                gene.setFieldValue("organism", organism);
                geneIdMap.put(identifier, gene);
            }
        }
        return gene;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected InterMineObject getMRNA(String mrnaIdentifier, String source, Organism organism, Model model) throws ObjectStoreException {
        // overriding getMRNA method to get more control over how and when the MRNA objects are stored
        InterMineObject mrna = null;
        if (mrnaIdMap.containsKey(mrnaIdentifier)) {
            // if mrnaIdMap contains the given mrnaIdentifier as a key then fetch the MRNA object
            // corresponding to the key
            mrna = mrnaIdMap.get(mrnaIdentifier);
        }
        else {
            // else create a new MRNA object and store in mrnaIdMap
            if (model.hasClassDescriptor(model.getPackageName() + ".MRNA")) {
                @SuppressWarnings("unchecked") Class<? extends InterMineObject> mrnaCls = (Class<? extends InterMineObject>) model.getClassDescriptorByName("MRNA").getType();
                mrna = getDirectDataLoader().createObject(mrnaCls);
                mrna.setFieldValue("primaryIdentifier", mrnaIdentifier);
                mrna.setFieldValue("source", source);
                mrna.setFieldValue("organism", organism);
                mrnaIdMap.put(mrnaIdentifier, mrna);
            }
        }
        return mrna;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void processFile(File file) {
        try {
            FileReader fileReader = new FileReader(file);
            BufferedReader reader = new BufferedReader(fileReader);

            System.err .println("reading " + getSequenceType() + " sequence from: " + file);

            SequenceIterator iter =
                    (SequenceIterator) SeqIOTools.fileToBiojava("fasta", getSequenceType(), reader);

            if (!iter.hasNext()) {
                System.err .println("no fasta sequences found - exiting");
                return;
            }

            while (iter.hasNext()) {
                Sequence bioJavaSequence = iter.nextSequence();
                processSequence(getOrganism(bioJavaSequence), bioJavaSequence);
            }
            // now storing all the Gene and MRNA objects
            storePendingObjectsToDataStore();
            reader.close();
            fileReader.close();
        } catch (BioException e) {
            throw new BuildException("sequence not in fasta format or wrong alphabet for: "
                    + file, e);
        } catch (NoSuchElementException e) {
            throw new BuildException("no fasta sequences in: " + file, e);
        } catch (FileNotFoundException e) {
            throw new BuildException("problem reading file - file not found: " + file, e);
        } catch (ObjectStoreException e) {
            throw new BuildException("ObjectStore problem while processing: " + file, e);
        } catch (IOException e) {
            throw new BuildException("error while closing FileReader for: " + file, e);
        }
    }

    /**
     * Stores all the created Gene and MRNA objects into the data store
     * @throws ObjectStoreException
     */
    private void storePendingObjectsToDataStore() throws ObjectStoreException {
        for (String geneIdentifier : geneIdMap.keySet()) {
            getDirectDataLoader().store(geneIdMap.get(geneIdentifier));
        }

        for (String mrnaIdentifier : mrnaIdMap.keySet()) {
            getDirectDataLoader().store(mrnaIdMap.get(mrnaIdentifier));
        }
    }
}

