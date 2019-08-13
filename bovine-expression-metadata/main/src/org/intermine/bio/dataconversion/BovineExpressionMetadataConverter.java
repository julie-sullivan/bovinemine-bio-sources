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
public class BovineExpressionMetadataConverter extends BioFileConverter
{
    //
    protected static final Logger LOG = Logger.getLogger(BovineExpressionMetadataConverter.class);

    private HashMap<String,Item> btoItems = new HashMap<String, Item>();
    private static final String DATASET_TITLE = "Metadata for Bovine RNASeq";
    private static final String DATA_SOURCE_NAME = "Metadata for Bovine RNASeq expression from USDA";
    private static final String TAXON_ID = "9913";
    private String orgRefId = getOrganism(TAXON_ID);

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public BovineExpressionMetadataConverter(ItemWriter writer, Model model) {
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
            if (Pattern.matches("Experiment", line[0])) {
                // skipping header
                continue;
            }

            String experiment = line[0];
            String tissue = line[1];
            String run = line[2];
            String biosample = line[3];
            String organSystem = line[4];
            String btoId = line[5];
            String btoName = line[6];
            String btoTopLevel = line[7];
            String btoAllLevels = line[8];
            String sex = line[9];
            String age = line[10];
            String breed = line[11];
            String individual = line[12];
            String lboId = line[13];
            String releaseDate = line[14];
            String spots = line[15];
            String bases = line[16];
            String spotsWithMates = line[17];
            String averageReadLength = line[18];
            String libraryName = line[19];
            String libraryStrategy = line[20];
            String librarySelection = line[21];
            String librarySource = line[22];
            String libraryLayout = line[23];
            String platform = line[24];
            String model = line[25];
            String sraStudy = line[26];
            String bioproject = line[27];
            String sraSample = line[28];
            String organismName = line[29];
            String sampleName = line[30];
            String sraSubmission = line[31];
            
            //String additionalTerms = line[37];
            Item item = createItem("ExpressionMetadata");
            if (!experiment.isEmpty()) {
                item.setAttribute("experiment", experiment);
            }
            else {
                System.out.println("experiment cannot be empty as it serves as a primaryIdentifier");
                System.exit(1);
            }

            item.setAttribute("experiment", experiment);
            item.setAttribute("tissue", tissue);
            item.setAttribute("run", run);
            item.setAttribute("biosample", biosample);
            item.setAttribute("organSystem", organSystem);
            item.setAttribute("btoId", btoId);
            item.setAttribute("btoName", btoName);
            item.setAttribute("btoTopLevel", btoTopLevel);
            item.setAttribute("btoAllLevels", btoAllLevels);
            item.setAttribute("sex", sex);
            item.setAttribute("age", age);
            item.setAttribute("breed", breed);
            item.setAttribute("individual", individual);
            item.setAttribute("lboId", lboId);
            item.setAttribute("releaseDate", releaseDate);
            item.setAttribute("spots", spots);
            item.setAttribute("bases", bases);
            item.setAttribute("spotsWithMates", spotsWithMates);
            item.setAttribute("averageReadLength", averageReadLength);
	    if (!libraryName.isEmpty()) {
                item.setAttribute("libraryName", libraryName);
            }
            //item.setAttribute("libraryName", libraryName);
            item.setAttribute("libraryStrategy", libraryStrategy);
            item.setAttribute("librarySelection", librarySelection);
            item.setAttribute("librarySource", librarySource);
            item.setAttribute("libraryLayout", libraryLayout);
            item.setAttribute("platform", platform);
            item.setAttribute("model", model);
            item.setAttribute("sraStudy", sraStudy);
            item.setAttribute("bioproject", bioproject);
            item.setAttribute("sraSample", sraSample);
            item.setAttribute("organismName", organismName);
            item.setAttribute("sampleName", sampleName);
            item.setAttribute("sraSubmission", sraSubmission);

            item.setReference("organism", getOrganism(TAXON_ID));

            if (!btoId.isEmpty()) {
                System.out.println("btoId: " + btoId);
                String[] btoIds = btoId.split(",");
                for (String btoIdentifier : btoIds) {
                    if (btoItems.containsKey(btoIdentifier)) {
                        item.addToCollection("brendaTissueOntology", btoItems.get(btoIdentifier).getIdentifier());
                    }
                    else {
                        Item btoItem = createItem("BRENDATerm");
                        btoItem.setAttribute("identifier", btoIdentifier);
                        item.addToCollection("brendaTissueOntology", btoItem.getIdentifier());
                        btoItem.addToCollection("samples", item.getIdentifier());
                        btoItems.put(btoIdentifier, btoItem);
                    }
                }
            }

            // Additional terms
//            if (!additionalTerms.isEmpty()) {
//                String[] brendaTissueOntologyTerms = additionalTerms.split(",");
//                for (String brendaTissueOntologyTerm : brendaTissueOntologyTerms) {
//                    String[] pair = brendaTissueOntologyTerm.split("\\|");
//                    if (btoItems.containsKey(pair[1])) {
//                        item.addToCollection("brendaTissueOntology", btoItems.get(pair[1]).getIdentifier());
//                    }
//                    else {
//                        Item btoItem = createItem("BRENDATerm");
//                        btoItem.setAttribute("identifier", pair[1]);
//                        btoItem.setAttribute("name", pair[0]);
//                        item.addToCollection("brendaTissueOntology", btoItem.getIdentifier());
//                        btoItem.addToCollection("samples", item.getIdentifier());
//                        btoItems.put(pair[1], btoItem);
//                    }
//                }
//            }

            try {
                store(item);
            } catch(Exception e) {
                System.out.println("Error while storing ExpressionMetadata item: " + item + "\nStacktrace: " + e);
            }
        }
    }

    /**
     *
     * {@inheritDoc}
     */
    @Override
    public void close() throws Exception {
        for (String key : btoItems.keySet()) {
            store(btoItems.get(key));
        }
    }
}
