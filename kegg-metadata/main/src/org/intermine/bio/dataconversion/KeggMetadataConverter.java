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
import java.io.Reader;
import java.lang.Exception;
import java.lang.String;
import java.lang.System;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;
import org.apache.log4j.Logger;
import org.intermine.util.FormattedTextParser;



/**
 * 
 * @author
 */
public class KeggMetadataConverter extends BioFileConverter
{
    //
    protected static final Logger LOG = Logger.getLogger(KeggMetadataConverter.class);
    private String orgRefId;
    private HashMap<String,Item> publicationItems = new HashMap<String, Item>();
    private HashMap<String,Item> goItems = new HashMap<String, Item>();
    private HashMap<String,Item> goAnnotation = new HashMap<String, Item>();
    private Map<String, String> goTerms = new HashMap<String, String>();
    private static final String DATASET_TITLE = "KEGG pathways metadata";
    private static final String DATA_SOURCE_NAME = "KEGG";

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public KeggMetadataConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {

     // assumes that the metadata file has unique entries
        Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);

        while (lineIter.hasNext()) {
            String[] line = lineIter.next();

         System.out.println(line.toString());

         String identifier = line[0];
            String pathwayClass = line[1];
            String pathwayMap = line[2];
            String koPathway = line[3];
             String pubmedId = line[4];
            String bsid = line[5];
            String description = line[6];
         //    String sampleTaxon = line[7];



           Item item = createItem("Pathway");
            if (!identifier.isEmpty()) {
                item.setAttribute("identifier", identifier);
               System.out.println(identifier);  
            }
            else {
                System.out.println("identifier cannot be empty as it serves as a primaryIdentifier");
                System.exit(1);
            }
            if (!pathwayClass.isEmpty()) { item.setAttribute("pathwayClass", pathwayClass); }
            if (!pathwayMap.isEmpty()) { item.setAttribute("pathwayMap", pathwayMap); }
            if (!koPathway.isEmpty()) { item.setAttribute("koPathway", koPathway); }
           if (!description.isEmpty()) { item.setAttribute("description", description); }
            if (!bsid.isEmpty()) {
                     String[] arraygo = bsid.split("--");
                      String bsid1 =arraygo[0];
                      String trimbsid = bsid1.substring(6);
			 item.setAttribute("bsid", trimbsid);  
           
                     String go =arraygo[1];
           if  (go.equals("PLACEHOLDER")){
           
            } else {                
                     String[] arr = go.split(" "); 
                     for (int i = 0; i < arr.length; i++) {
                          String goterm = arr[i];
                       System.out.println(goterm);                     
            if (!goterm.isEmpty()) {
            
                 String goTermRefId = getGoTerm(goterm.trim());

                 Item goAnnotation = createItem("GOAnnotation");
                  goAnnotation.setReference("ontologyTerm", goTermRefId);

                  goAnnotation.setReference("subject", item);
                  item.addToCollection("goAnnotation", goAnnotation);
   

            try {
                store(goAnnotation);
            } catch (Exception e) {
                System.out.println("Error while storing goAnnotation item: " + item + "\nStacktrace: " + e);
            }

}
 }              
}  
}
         //   String TAXON_ID = sampleTaxon;
         //   orgRefId = getOrganism(TAXON_ID);
        //    item.setReference("organism", getOrganism(TAXON_ID));

           if (!pubmedId.isEmpty()) {
                   String[] array = pubmedId.split(" , ");
                  for (int i = 0; i < array.length; i++) {
                          String pubMed = array[i];
                if (publicationItems.containsKey(pubMed)) {
                    item.addToCollection("publications", publicationItems.get(pubMed).getIdentifier());
                }
                else {
                    Item pubItem = createItem("Publication");
                    pubItem.setAttribute("pubMedId", pubMed);
                   item.addToCollection("publications", pubItem.getIdentifier());
                    publicationItems.put(pubMed, pubItem);
                }
              }
            }

             try {
                store(item);	
            } catch(Exception e) {
                System.out.println("Error while storing ExpressionMetadata item: " + item + "\nStacktrace: " + e);
            }

       }
    }

   @Override
    public void close() throws Exception {
        for (String key : publicationItems.keySet()) {
            store(publicationItems.get(key));
        }
    }

    private String getGoTerm(String identifier) {
        String refId = goTerms.get(identifier);
        if (refId == null) {
            Item item = createItem("GOTerm");
            item.setAttribute("identifier", identifier);
            refId = item.getIdentifier();
            goTerms.put(identifier, refId);
            try {
                store(item);
            } catch (Exception e) {
               System.out.println("Error while storing GOTerm item: " + item + "\nStacktrace: " + e);
            }
        }
        return refId;
    }



}






