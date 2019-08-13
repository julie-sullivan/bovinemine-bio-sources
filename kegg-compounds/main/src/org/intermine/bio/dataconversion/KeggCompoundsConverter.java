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
public class KeggCompoundsConverter extends BioFileConverter
{
    //
    protected static final Logger LOG = Logger.getLogger(KeggCompoundsConverter.class);
    private String orgRefId;
    private HashMap<String,Item> compoundsItems = new HashMap<String, Item>();
    private HashMap<String,Item> compsItems = new HashMap<String, Item>();
    private static final String DATASET_TITLE = "KEGG pathways metadata";
    private static final String DATA_SOURCE_NAME = "KEGG";

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public KeggCompoundsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {

        Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);

        while (lineIter.hasNext()) {
            String[] line = lineIter.next();

         System.out.println(line.toString());

         String identifier = line[0];
            String compounds = line[1];
         



           Item item = createItem("Pathway");
            if (!identifier.isEmpty()) {
                item.setAttribute("identifier", identifier);
               System.out.println(identifier);  
            }
            else {
                System.out.println("identifier cannot be empty as it serves as a primaryIdentifier");
                System.exit(1);
            }
              if (!compounds.isEmpty()) {
                   String[] array = compounds.split(" ");
                  for (int i = 0; i < array.length; i++) {
                          String onecomp = array[i];
                if (compoundsItems.containsKey(onecomp)) {
                    item.addToCollection("compound", compoundsItems.get(onecomp).getIdentifier());
                }
                else {
                    Item compsItems = createItem("Compound");
                    compsItems.setAttribute("primaryIdentifier", onecomp);
                   item.addToCollection("compound", compsItems.getIdentifier());
                    compoundsItems.put(onecomp, compsItems);
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
        for (String key : compoundsItems.keySet()) {
            store(compoundsItems.get(key));
        }
    }

}






