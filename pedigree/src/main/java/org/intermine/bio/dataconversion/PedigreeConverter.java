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
public class PedigreeConverter extends BioFileConverter
{
    //
    private static final Logger LOG = Logger.getLogger(PedigreeConverter.class);
    private static final String PROP_FILE = "config.properties";
    private Set<String> taxonIds;
    private static final String DATASET_TITLE = "Pedigree data set";
    private static final String DATA_SOURCE_NAME = "UNKNOWN";
    private Map<String, String> genes = new HashMap<String, String>();
    private Set<String> relates = new HashSet<String>();
     protected IdResolver rslv = null;
    private Map<String, String> configs = new HashMap<String, String>();
  private static String evidenceRefId = null;
   

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public PedigreeConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
       readConfig();
    }

//
     public void setPedigreeOrganisms(String taxonIds) {
        this.taxonIds = new HashSet<String>(Arrays.asList(StringUtils.split(taxonIds, " ")));
    }

    /**
     * Sets the list of organisms to process only if the genes are homologues for organism of
     * interest.  Otherwise ignore.
     * @param taxonIds list of taxon IDs to process
     */
    public void setPedigreeRelates(String taxonIds) {
        this.relates = new HashSet<String>(Arrays.asList(StringUtils.split(taxonIds, " ")));
    }

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
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {

        // init resolver
    
        if (taxonIds == null || taxonIds.isEmpty()) {
            throw new IllegalArgumentException("No organism data provided for Ensembl Compara");
        }
        File file = getCurrentFile();
        if (file == null) {
            throw new FileNotFoundException("No valid data files found.");
        }
        String fileName = file.getName();
        String[] bits = fileName.split("_");
        boolean processFile = false;
        for (String bit : bits) {
            if (taxonIds.contains(bit)) {
                processFile = true;
            } else if (!relates.isEmpty() && !relates.contains(bit)) {
                // this file contains an organism not listed in the project XML file
                return;
            }
        }
        if (!processFile) {
            return;
        }

        String lastGene1 = "";
        String lastGene2 = "";
        Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length < 3 && StringUtils.isNotEmpty(line.toString())) {
                throw new RuntimeException("Invalid line, should be 2 columns but is '"
                        + line.length + "' instead");
            }

            String gene1 = line[1];
            String gene2 = line[2];
            String family = line[0];
            String relate = line[3];
            String relate1 = line[4];
         
            

            if (gene1.equals(lastGene1) && gene2.equals(lastGene2)) {
                // file isn't unique
                continue;
            }

            String refId2 = parseIndividual(bits[1], gene2);
            String refId1 = parseIndividual(bits[0], gene1);
            if (refId1 == null || refId2 == null) {
                continue;
            }
            processRelates(refId1, refId2, family, relate);
            processRelates2(refId2, refId1, family, relate1);
            lastGene1 = gene1;
            lastGene2 = gene2;
        }
    }


       private void processRelates(String gene1, String gene2, String family, String relate)
        throws ObjectStoreException {
        Item rel = createItem("Relationships");
        rel.setReference("relationships", gene1);
        rel.setReference("individual", gene2);
        rel.setReference("relationshiplookup", getRelate(relate));
        String refId = rel.getIdentifier();
        System.out.println(">>>>>>>>>>>>>>"+ refId);
        rel.setAttribute("type", family);
        store(rel);
    }
    private void processRelates2(String gene1, String gene2, String family, String relate1)
        throws ObjectStoreException {
        Item rel = createItem("Relationships");
        rel.setReference("relationships", gene1);
        rel.setReference("individual", gene2);
        rel.setReference("relationshiplookup", getRelate(relate1));
        String refId = rel.getIdentifier();
        System.out.println(">>>>>>>>>>>>>>"+ refId);
        rel.setAttribute("type", family);
        store(rel);
    }






    private String parseIndividual(String taxonId, String identifier)
        throws ObjectStoreException {
        if (StringUtils.isBlank(identifier)) {
            return null;
        }
        String newIdentifier = identifier;
       String refId = genes.get(newIdentifier);
        if (refId == null) {
            String fieldName = getConfig(taxonId);
            if (fieldName == null) {
                throw new IllegalArgumentException("no config found");
            }
            Item item = createItem("Individual");
            item.setAttribute(fieldName, newIdentifier);
            item.setReference("organism", getOrganism(taxonId));
            store(item);
            refId = item.getIdentifier();
            genes.put(newIdentifier, refId);
        }
        return refId;
    }

    private String getConfig(String taxonId) {
        String identifierField = configs.get(taxonId);
        if (identifierField == null) {
            identifierField  = configs.get("default");
        }
        return identifierField;
    }
  
    

    private String getRelate(String identifier) {
          String newId = identifier;
         Item relate = createItem("Relationshiplookup"); 
          relate.setAttribute("primaryIdentifier", newId);     
         try{ 
         store(relate);
         } catch(Exception e) {
         System.out.println("Error while storing relate item: " + relate + "\nStacktrace: " + e);
         }
         String refId2 = relate.getIdentifier();
         return refId2;
       
    }





    /**
     * 
     *
     * {@inheritDoc}
     */
}
