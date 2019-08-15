package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2016 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Reader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;


/**
 * 
 * @author
 */
public class TreefamExtendedConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(TreefamExtendedConverter.class);
    private static final String PROP_FILE = "treefam-extended_config.properties";
    private static final String DATASET_TITLE = "TreeFam data set";
    private static final String DATA_SOURCE_NAME = "TreeFam";
    private static final String EVIDENCE_CODE_ABBR = "AA";
    private static final String EVIDENCE_CODE_NAME = "Amino acid sequence comparison";
    private HashMap<String, String> organismMap = new HashMap<String, String>();

    private Set<String> homologues = new HashSet<String>();
    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> configs = new HashMap<String, String>();
    private static String evidenceRefId = null;
    protected String version = null;
    private String geneSource = null;

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public TreefamExtendedConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        readConfig();
        organismMap.put("Bos taurus", "9913");
        organismMap.put("Canis familiaris", "9615");
        organismMap.put("Homo sapiens", "9606");
        organismMap.put("Mus musculus", "10090");
        organismMap.put("Rattus norvegicus", "10116");
        organismMap.put("Equus caballus", "9796");
        organismMap.put("Gallus gallus", "9031");
        organismMap.put("Sus scrofa", "9823");
    }

    /**
     * Set the version for TreeFam. This value is used in 'source' attribute
     * of Gene for integration
     * @param version the version
     */
    public void setVersion(String version) {
        this.version = version;
    }

    public void setGeneSource(String geneSource) {
        this.geneSource = geneSource;
    }

    public String getGeneSource() {
        return this.geneSource;
    }

    /**
     * Read config file
     */
    private void readConfig() {
        Properties props = new Properties();
        try {
            props.load(getClass().getClassLoader().getResourceAsStream(PROP_FILE));
        } catch (IOException e) {
            throw new RuntimeException("Problem loading properties '" + PROP_FILE + "'", e);
        }

        for (Map.Entry<Object, Object> entry: props.entrySet()) {
            String key = (String) entry.getKey();
            String value = ((String) entry.getValue()).trim();
            if (StringUtils.isEmpty(key) || StringUtils.isEmpty(value)) {
                throw new RuntimeException("Problem loading properties '" + PROP_FILE + "' on line "
                        + key + " = " + value);
            }
            configs.put(key, value);
        }
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {

        if (version == null) {
            throw new IllegalArgumentException("No version provided for Ensembl");
        }

        String lastGene1 = "";
        String lastGene2 = "";
        Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length < 7 && StringUtils.isNotEmpty(line.toString())) {
                throw new RuntimeException("Invalid line, should be 7 columns but is " + line.length + " instead");
            }

            String gene1 = line[0];
            String organism1 = line[1];
            String gene2 = line[2];
            String organism2 = line[3];
            String lca = line[4];
            String homologyType = line[5];
            // TODO: what is line[6]?

            if (gene1.equals(lastGene1) && gene2.equals(lastGene2)) {
                // file isn't unique
                continue;
            }

            String refId1 = parseGene(gene1, organism1);
            String refId2 = parseGene(gene2, organism2);

            if (refId1 == null || refId2 == null) {
                continue;
            }

            processHomologue(refId1, refId2, lca, homologyType);
            processHomologue(refId2, refId1, lca, homologyType);
            lastGene1 = gene1;
            lastGene2 = gene2;
        }
    }

    /**
     * Creates a Gene item for a given identifier and stores it in the database
     * @param identifier
     * @param organismName
     * @return
     * @throws ObjectStoreException
     */
    private String parseGene(String identifier, String organismName) throws ObjectStoreException {
        if (StringUtils.isBlank(identifier)) {
            return null;
        }
        String refId = genes.get(identifier);
        if (refId == null) {
            String taxonId = organismMap.get(organismName);
            String fieldName = getConfig(taxonId);
            if (fieldName == null) {
                throw new IllegalArgumentException("no config found");
            }
            Item item = createItem("Gene");
            item.setAttribute(fieldName, identifier);
	    /*
            if ("9913".equals(taxonId) || "9940".equals(taxonId)) {
                item.setAttribute("source", version);
            }
	    */
	    item.setAttribute("source", this.geneSource);
            item.setReference("organism", getOrganism(taxonId));
            store(item);
            refId = item.getIdentifier();
            genes.put(identifier, refId);
        }
        return refId;
    }

    /**
     * Creates a Homologue item that represents the homology between two Gene items
     * @param gene1
     * @param gene2
     * @param lastCommonAncestor
     * @param homologyType
     * @throws ObjectStoreException
     */
    private void processHomologue(String gene1, String gene2, String lastCommonAncestor, String homologyType) throws ObjectStoreException {
        Item homologue = createItem("Homologue");
        homologue.setReference("gene", gene1);
        homologue.setReference("homologue", gene2);
        homologue.addToCollection("evidence", getEvidence());
        homologue.setAttribute("type", homologyType.contains("ortholog") ? "orthologue" : "paralogue");
        homologue.setAttribute("lastCommonAncestor", lastCommonAncestor);
        homologue.setAttribute("homologyType", homologyType);
        store(homologue);
    }

    /**
     * Get configuration from config for a given taxonId
     * @param taxonId
     * @return
     */
    private String getConfig(String taxonId) {
        String identifierField = configs.get(taxonId);
        if (identifierField == null) {
            identifierField  = configs.get("default");
        }
        return identifierField;
    }

    /**
     * Get Evidence item for the current data source
     * @return
     * @throws ObjectStoreException
     */
    private String getEvidence() throws ObjectStoreException {
        if (evidenceRefId == null) {
            Item item = createItem("OrthologueEvidenceCode");
            item.setAttribute("abbreviation", EVIDENCE_CODE_ABBR);
            item.setAttribute("name", EVIDENCE_CODE_NAME);
            try {
                store(item);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }
            String refId = item.getIdentifier();

            item = createItem("OrthologueEvidence");
            item.setReference("evidenceCode", refId);
            try {
                store(item);
            } catch (ObjectStoreException e) {
                throw new ObjectStoreException(e);
            }

            evidenceRefId = item.getIdentifier();
        }
        return evidenceRefId;
    }

}
