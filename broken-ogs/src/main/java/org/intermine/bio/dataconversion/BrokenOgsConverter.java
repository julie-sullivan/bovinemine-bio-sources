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

import java.io.Reader;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;
import org.intermine.util.FormattedTextParser;

import java.lang.NumberFormatException;
import java.lang.System;
import java.util.Iterator;
import java.util.ArrayList;

/**
 * 
 * @author
 */
public class BrokenOgsConverter extends BioFileConverter
{
    //
    private static final String DATASET_TITLE = "btau_OGSv2";
    private static final String DATA_SOURCE_NAME = "Bos taurus Official Gene Set OGSv2 (on UMD3.1)";
    private static final String TAXON_ID = "9913";
    private static final String SOURCE = "btau_OGSv2";
    private String orgRefId;
    private ArrayList<Item> itemsToStore = new ArrayList<Item>();
    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public BrokenOgsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        orgRefId = getOrganism(TAXON_ID);
        Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            System.out.println(line[0] + "\t" + line[1]);
            String geneIdentifier = line[0].trim();
            String[] transcriptIdentifierList = line[1].split(",");
            Item geneItem = createItem("Gene");
            geneItem.setAttribute("primaryIdentifier", geneIdentifier);
            geneItem.setAttribute("source", SOURCE);
            geneItem.setAttribute("status", "failed liftOver");
            geneItem.setReference("organism", orgRefId);
            itemsToStore.add(geneItem);
            for (String transcriptIdentifier : transcriptIdentifierList) {
                Item transcriptItem = createItem("MRNA");
                transcriptItem.setAttribute("primaryIdentifier", transcriptIdentifier);
                transcriptItem.setAttribute("status", "failed liftOver");
                transcriptItem.setReference("organism", orgRefId);
                transcriptItem.setReference("gene", geneItem.getIdentifier());
                geneItem.addToCollection("transcripts", transcriptItem);
                itemsToStore.add(transcriptItem);
            }

            for (Item item : itemsToStore) {
                try {
                    store(item);
                } catch (Exception e) {
                    System.out.println("Error while storing item: " + item + "\n" + "Exception stack trace: " + e);
                }
            }
            itemsToStore.clear();
        }
    }
}
