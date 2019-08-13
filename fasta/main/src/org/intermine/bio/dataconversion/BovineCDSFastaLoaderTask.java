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

import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.tools.ant.BuildException;

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
 * A fasta loader that understand the headers of RefSeq CDS fasta files and can make the
 * appropriate extra objects and references.
 * @author
 */
public class BovineCDSFastaLoaderTask extends BovineFeatureFastaLoaderTask
{

    protected static final String HEADER_REGEX = ".(\\S+)\\s(\\S+)";
    private String classAttribute = "primaryIdentifier";
    private String suffix = "-CDS";
    private int storeCount = 0;

    /**
     * {@inheritDoc}
     */
    @Override
    protected void processSequence(Organism organism, Sequence bioJavaSequence) throws ObjectStoreException {

        if (organism == null) {
            return;
        }

        org.intermine.model.bio.Sequence proteinSequenceObject = getDirectDataLoader().createObject(org.intermine.model.bio.Sequence.class);

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
        String attributeValue;
        if (sequenceName.endsWith(suffix)) {
            attributeValue = sequenceName;
        }
        else {
            attributeValue = sequenceName + suffix;
        }

        // set primaryIdentifier
        try {
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
        Annotation annotation = bioJavaSequence.getAnnotation();
        String mrnaIdentifier = bioJavaSequence.getName();
        String header = (String) annotation.getProperty("description");
        Pattern p = Pattern.compile(HEADER_REGEX);
        Matcher m = p.matcher(header);
        String proteinIdentifier = "";
        if (m.matches()) {
            proteinIdentifier = m.group(2);
        }
        if (proteinIdentifier != "") {
            bioEntity.setFieldValue("proteinIdentifier", proteinIdentifier);
        }
        ObjectStore os = getIntegrationWriter().getObjectStore();
        Model model = os.getModel();
        if (model.hasClassDescriptor(model.getPackageName() + ".CDS")) {
            Class<? extends FastPathObject> cdsCls =
            model.getClassDescriptorByName("CDS").getType();
            if (!DynamicUtil.isInstance(bioEntity, cdsCls)) {
                throw new RuntimeException("the InterMineObject passed to "
                + "BovineCDSFastaLoaderTask.extraProcessing() is not a "
                + "CDS: " + bioEntity);
            }
            InterMineObject mrna = getMRNA(mrnaIdentifier, getGeneSource(), organism, model);
            if (mrna != null) {
                bioEntity.setFieldValue("transcript", mrna);
                mrna.setFieldValue("proteinIdentifier",proteinIdentifier);
               // bioEntity.setFieldValue("mrnaIdentifier", mrnaIdentifier);
               // bioEntity.setFieldValue("mrna", mrna);
		//try {
                  //  HashSet polypeptideCollection = (HashSet) mrna.getFieldValue("polypeptide");
                    //polypeptideCollection.add(bioEntity);
                   // mrna.setFieldValue("polypeptide", polypeptideCollection);
                   // mrnaIdMap.put(mrnaIdentifier, mrna);
               // } catch(IllegalAccessException e) {
                 //   e.printStackTrace();
               // }
            }
            //Location loc = getLocationFromHeader(header, (SequenceFeature) bioEntity, organism);
            //getDirectDataLoader().store(loc);
        } else {
            throw new RuntimeException("Trying to load CDS sequence but CDS does not exist in the data model");
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected String getIdentifier(Sequence bioJavaSequence) {
        Annotation annotation = bioJavaSequence.getAnnotation();
        String mrnaIdentifier = bioJavaSequence.getName();
        return mrnaIdentifier;
    }
}
