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
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;
import org.intermine.util.FormattedTextParser;

import java.util.Map;
import java.util.HashMap;
import java.lang.NumberFormatException;
import java.lang.System;
import java.util.Iterator;


/**
 * A Converter that loads cross references, across gene sets, for Gene entity.
 * @author
 */
public class CrossReferencesConverter extends BioFileConverter
{

    private static final String DATASET_TITLE = null;
    private static final String DATA_SOURCE_NAME = null;
    Map<String, Item> geneItemMap = new HashMap<String, Item>();
    String taxonId;
    String organismReferenceId;

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public CrossReferencesConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * setter for taxonId, when its provided by project.xml
     * @param taxonId
     */
    public void setTaxonId(String taxonId) {
        this.taxonId = taxonId;
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        File currentFile = getCurrentFile();
        organismReferenceId = getOrganism(taxonId);
        Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);

        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            String[] referrerInfo = line[0].trim().split(":");
            String referrerPrimaryIdentifier = referrerInfo[0];
            String referrerSource = referrerInfo[1];

            Item geneItem1;
            if (geneItemMap.containsKey(referrerPrimaryIdentifier)) {
                geneItem1 = geneItemMap.get(referrerPrimaryIdentifier);
            }
            else {
                geneItem1 = createItem("Gene");
                geneItem1.setAttribute("primaryIdentifier", referrerPrimaryIdentifier);
                geneItem1.setAttribute("source", referrerSource);
                geneItem1.setReference("organism", organismReferenceId);
                geneItemMap.put(referrerPrimaryIdentifier, geneItem1);
            }

            String[] refereeEntryList = line[1].trim().split(",");
            for (String refereeEntry : refereeEntryList) {
                String[] refereeInfo = refereeEntry.split(":");
                String refereePrimaryIdentifier = refereeInfo[0];
                String refereeSource = refereeInfo[1];
                System.out.println("REFERRER: " + referrerPrimaryIdentifier + " <> REFERREE: " + refereePrimaryIdentifier);

                Item geneItem2;
                if (geneItemMap.containsKey(refereePrimaryIdentifier)) {
                    geneItem2 = geneItemMap.get(refereePrimaryIdentifier);
                }
                else {
                    geneItem2 = createItem("Gene");
                    geneItem2.setAttribute("primaryIdentifier", refereePrimaryIdentifier);
                    geneItem2.setAttribute("source", refereeSource);
                    geneItem2.setReference("organism", organismReferenceId);
                    geneItemMap.put(refereePrimaryIdentifier, geneItem2);
                }

                // create xRef linking gene1 and gene2
                Item xrefItem = createItem("xRef");
                xrefItem.setAttribute("refereeSource", refereeSource);
                xrefItem.setReference("organism", organismReferenceId);
                xrefItem.setReference("referrer", geneItem1.getIdentifier());
                xrefItem.setReference("referee", geneItem2.getIdentifier());

                           
                // add xref to dbCrossReferences collection of gene1
                String xrefDbIdentifier = xrefItem.getIdentifier();
                geneItem1.addToCollection("dbCrossReferences", xrefDbIdentifier);
                geneItemMap.put(referrerPrimaryIdentifier, geneItem1);
                storeItem(xrefItem);
            }
        }
    }

    /**
     * Store a given item
     * @param item
     */
    protected void storeItem(Item item) {
        try {
            store(item);
        } catch (Exception e) {
            System.out.println("Error while storing item: " + item);
            System.out.println("Exception stacktrace: " + e);
        }
    }

    /**
     * Store all items in a given Map
     * @param itemMap
     */
    protected void storeAllItems(Map<String, Item> itemMap) {
        for (String key : itemMap.keySet()) {
            storeItem(itemMap.get(key));
        }
    }

    /**
     *
     * {@inheritDoc}
     */
    @Override
    public void close() throws Exception {
        storeAllItems(geneItemMap);
    }
}
