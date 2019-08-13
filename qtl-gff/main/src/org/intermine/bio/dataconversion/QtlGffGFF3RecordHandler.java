package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2015 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */


import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.metadata.Model;
import org.intermine.metadata.StringUtil;
import org.intermine.xml.full.Item;

import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.intermine.bio.util.OrganismRepository;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.SAXParser;
import org.intermine.xml.full.ReferenceList;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;


/**
 * A converter/retriever for the QtlGff dataset via GFF files.
 */

public class QtlGffGFF3RecordHandler extends GFF3RecordHandler
{

    private static final Logger LOG = Logger.getLogger(QtlGffGFF3RecordHandler.class);
    private static final String CLINICAL_MEASUREMENT_ONTOLOGY = "Clinical Measurement Ontology";
    private static final String VERTEBRATE_TRAIT_ONTOLOGY = "Vertebrate Trait Ontology";
    private static final String LIVESTOCK_PRODUCT_TRAIT_ONTOLOGY = "Livestock Product Trait Ontology";
    private String geneSource = null;
    HashMap<String, Item> delayedItems = new HashMap<String, Item>();

    /**
     * Create a new QtlGffGFF3RecordHandler for the given data model.
     * @param model the model for which items will be created
     */
    public QtlGffGFF3RecordHandler (Model model) {
        super(model);
    }

    public void setGeneSource(String geneSource) {
        this.geneSource = geneSource;
    }

    public String getGeneSource() {
        return this.geneSource;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(GFF3Record record) { 
        Item feature = getFeature();
        String clsName = feature.getClassName();
	//System.out.println(record.getAttributes());
        if(clsName.equals("QTL")) {
            if (record.getAttributes().get("ID") != null) {
                String primaryIdentifier = record.getAttributes().get("ID").iterator().next();
                feature.setAttribute("primaryIdentifier", primaryIdentifier);
            }
            if (record.getAttributes().get("QTL_ID") != null) {
                String qtlId = record.getAttributes().get("QTL_ID").iterator().next();
                feature.setAttribute("qtlId", qtlId);
            }
            if (record.getAttributes().get("Name") != null) {
                String name = record.getAttributes().get("Name").iterator().next();
                feature.setAttribute("name", name);
            }
            if (record.getAttributes().get("qtl_type") != null) {
                String qtl_type = record.getAttributes().get("qtl_type").iterator().next().replace("_"," ");
                feature.setAttribute("type", qtl_type);
            }
            if (record.getAttributes().get("trait_ID") != null) {
                String traitName = record.getAttributes().get("trait_ID").iterator().next();
                feature.setAttribute("traitId", traitName);
            }
            if (record.getAttributes().get("trait") != null) {
                String traitName = record.getAttributes().get("trait").iterator().next();
                feature.setAttribute("trait", traitName);
            }
            if (record.getAttributes().get("Abbrev") != null) {
                String abbrevName = record.getAttributes().get("Abbrev").iterator().next().replace("-","%");
                feature.setAttribute("abbreviation", abbrevName);
            }
            if (record.getAttributes().get("Additive_Effect") != null) {
                String addEffect = record.getAttributes().get("Additive_Effect").iterator().next();
                feature.setAttribute("additiveEffect", addEffect);
            }
            if (record.getAttributes().get("Bayes-value") != null) {
                String bayesValue = record.getAttributes().get("Bayes-value").iterator().next();
                feature.setAttribute("bayesValue", bayesValue);
            }
            if (record.getAttributes().get("Dominance_Effect") != null) {
                String dominanceEffect = record.getAttributes().get("Dominance_Effect").iterator().next();
                feature.setAttribute("dominanceEffect", dominanceEffect);
            }
            if (record.getAttributes().get("breed") != null) {
                String breed = record.getAttributes().get("breed").iterator().next();
                feature.setAttribute("breed", breed);
            }
            if (record.getAttributes().get("Map_Type") != null) {
                String mapType = record.getAttributes().get("Map_Type").iterator().next();
                feature.setAttribute("mapType", mapType);
            }
            if (record.getAttributes().get("Model") != null) {
                String model = record.getAttributes().get("Model").iterator().next();
                feature.setAttribute("model", model);
            }
            if (record.getAttributes().get("Test_Base") != null) {
                String testBase = record.getAttributes().get("Test_Base").iterator().next();
                feature.setAttribute("testBase", testBase);
            }
            if (record.getAttributes().get("Significance") != null) {
                String significanceValue = record.getAttributes().get("Significance").iterator().next();
                feature.setAttribute("significance", significanceValue);
            }
            if (record.getAttributes().get("P-value") != null) {
                String pValue = record.getAttributes().get("P-value").iterator().next();
                feature.setAttribute("pValue", pValue);
            }
            //if (record.getAttributes().get("Variance") != null) {
                //String varianceValue = record.getAttributes().get("Variance").iterator().next();
                //feature.setAttribute("variance", varianceValue);
           //}
            if (record.getAttributes().get("peak_cM") != null) {
                String peakCmValue = record.getAttributes().get("peak_cM").iterator().next();
                feature.setAttribute("peakCentimorgan", peakCmValue);
            }
            if (record.getAttributes().get("Likelihood_Ratio") != null) {
                String likelihoodRatioValue = record.getAttributes().get("Likelihood_Ratio").iterator().next();
                feature.setAttribute("likelihoodRatio", likelihoodRatioValue);
            }
            if (record.getAttributes().get("LOD-score") != null) {
                String lodScoreValue = record.getAttributes().get("LOD-score").iterator().next();
                feature.setAttribute("lodScore", lodScoreValue);
            }
            if (record.getAttributes().get("LS-means") != null) {
                String lsMeansValue = record.getAttributes().get("LS-means").iterator().next();
                feature.setAttribute("lsMeans", lsMeansValue);
            }
            if (record.getAttributes().get("CMO_name") != null) {
                String cmoName = record.getAttributes().get("CMO_name").iterator().next();
                Item cmoTerm = getOntologyTerm(cmoName, CLINICAL_MEASUREMENT_ONTOLOGY);
                feature.setReference("cmoName", cmoTerm.getIdentifier());
            }
            if (record.getAttributes().get("VTO_name") != null) {
                String vtoName = record.getAttributes().get("VTO_name").iterator().next();
                Item vtoTerm = getOntologyTerm(vtoName, VERTEBRATE_TRAIT_ONTOLOGY);
                feature.setReference("vtoName", vtoTerm.getIdentifier());
            }
            if (record.getAttributes().get("PTO_name") != null) {
                String ptoName = record.getAttributes().get("PTO_name").iterator().next();
                Item ptoTerm = getOntologyTerm(ptoName, LIVESTOCK_PRODUCT_TRAIT_ONTOLOGY);
                feature.setReference("ptoName", ptoTerm.getIdentifier());
            }
            if (record.getAttributes().get("PUBMED_ID") != null) {
                String pubMedId = record.getAttributes().get("PUBMED_ID").iterator().next();
                Item publication = getPublication(pubMedId);
                feature.setReference("publication", publication.getIdentifier());
            }
            // more attributes
            if (record.getAttributes().get("gene_ID") != null) {
                String geneIdentifier = record.getAttributes().get("gene_ID").iterator().next();
                // assuming that all genes are NCBI RefSeq genes
                Item gene = getGene(geneIdentifier, "RefSeq");
                //Item gene = getGene(geneIdentifier, this.geneSource);
                feature.setReference("gene", gene.getIdentifier());
            }
            if (record.getAttributes().get("FlankMarkers") != null) {
                List<String> markers = record.getAttributes().get("FlankMarkers");
                feature.setAttribute("flankMarkers", StringUtils.join(markers, ", "));
                ArrayList<String> sequenceAlterations = new ArrayList<String>();

                for (String marker : markers) {
                    if (marker.startsWith("rs")) {
                        Item sequenceAlteration = getSequenceAlteration(marker);
                        sequenceAlterations.add(sequenceAlteration.getIdentifier());
                    }
                }
                for (int i = 0; i < sequenceAlterations.size(); i++) {
                    feature.addToCollection("snpAsFlankMarkers", sequenceAlterations.get(i));
                }
            }
        }
    }

    /**
     * Get an Item representation of a Gene
     * @param identifier
     * @param source
     * @return
     */
    private Item getGene(String identifier, String source) {
        Item gene = null;
        if (delayedItems.containsKey(identifier)) {
            gene = delayedItems.get(identifier);
        }
        else {
            gene = converter.createItem("Gene");
            gene.setAttribute("primaryIdentifier", identifier);
            gene.setAttribute("source", source);
            delayedItems.put(identifier, gene);
        }
        return gene;
    }

    /**
     * Get an Item representation of a SequenceAlteration
     * @param identifier
     * @return
     */
    private Item getSequenceAlteration(String identifier) {
        Item sequenceAlteration = null;
        if (delayedItems.containsKey(identifier)) {
            sequenceAlteration = delayedItems.get(identifier);
        }
        else {
            sequenceAlteration = converter.createItem("SequenceAlteration");
            sequenceAlteration.setAttribute("primaryIdentifier", identifier);
            delayedItems.put(identifier, sequenceAlteration);
        }
        return sequenceAlteration;
    }

    /**
     * Get an Item representation of a subclass of OntologyTerm based on ontologyName
     * @param termName
     * @param ontologyName
     * @return
     */
    private Item getOntologyTerm(String termName, String ontologyName) {
        Item ontologyTerm = null;
        String key = ontologyName + ":" + termName;
        if (delayedItems.containsKey(key)) {
            ontologyTerm = delayedItems.get(key);
        }
        else {
            if (ontologyName.equals(CLINICAL_MEASUREMENT_ONTOLOGY)) {
                ontologyTerm = converter.createItem("CMOTerm");
                ontologyTerm.setAttribute("name", termName);
                Item ontology = getOntology(ontologyName);
                ontologyTerm.setReference("ontology", ontology);
            }
            else if (ontologyName.equals(VERTEBRATE_TRAIT_ONTOLOGY)) {
                ontologyTerm = converter.createItem("VTTerm");
                ontologyTerm.setAttribute("name", termName);
                Item ontology = getOntology(ontologyName);
                ontologyTerm.setReference("ontology", ontology);
            }
            else if (ontologyName.equals(LIVESTOCK_PRODUCT_TRAIT_ONTOLOGY)) {
                ontologyTerm = converter.createItem("LPTTerm");
                ontologyTerm.setAttribute("name", termName);
                Item ontology = getOntology(ontologyName);
                ontologyTerm.setReference("ontology", ontology);
            }
            delayedItems.put(key, ontologyTerm);
        }
        return ontologyTerm;
    }

    /**
     * Get an Item representation of an Ontology
     * @param ontologyName
     * @return
     */
    private Item getOntology(String ontologyName) {
        Item ontology = null;
        if (delayedItems.containsKey(ontologyName)) {
            ontology = delayedItems.get(ontologyName);
        }
        else {
            ontology = converter.createItem("Ontology");
            ontology.setAttribute("name", ontologyName);
            if (ontologyName.equals(CLINICAL_MEASUREMENT_ONTOLOGY)) {
                ontology.setAttribute("url", "https://bioportal.bioontology.org/ontologies/CMO");
            }
            else if (ontologyName.equals(VERTEBRATE_TRAIT_ONTOLOGY)) {
                ontology.setAttribute("url", "https://bioportal.bioontology.org/ontologies/VT");
            }
            else if (ontologyName.equals(LIVESTOCK_PRODUCT_TRAIT_ONTOLOGY)) {
                ontology.setAttribute("url", "https://bioportal.bioontology.org/ontologies/LPT");
            }
            delayedItems.put(ontologyName, ontology);
        }
        return ontology;
    }
  
    /**
     * Get an Item representation of a Publication
     * @param pubMedId
     * @return
     */
    private Item getPublication(String pubMedId) {
        Item publication = null;
        if (delayedItems.containsKey(pubMedId)) {
            publication = delayedItems.get(pubMedId);
        }
        else {
            publication = converter.createItem("Publication");
            publication.setAttribute("pubMedId", pubMedId);
            delayedItems.put(pubMedId, publication);
        }
        return publication;
    }

    /**
     * Returns all delayed items
     * @return
     */
    @Override
    public Collection<Item> getFinalItems() {
        ArrayList<Item> finalItems = new ArrayList<Item>();
        for (String key : delayedItems.keySet()) {
            finalItems.add(delayedItems.get(key));
        }
        return finalItems;
    }
}
