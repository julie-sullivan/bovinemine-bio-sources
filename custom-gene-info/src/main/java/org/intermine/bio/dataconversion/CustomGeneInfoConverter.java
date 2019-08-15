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

import java.lang.NumberFormatException;
import java.lang.System;
import java.util.Iterator;

/**
 * A Converter that loads NCBI Gene ID to Ensembl Gene ID mapping to enable
 * proper integration of datasets using heterogenous identifiers.
 *
 * @author Deepak Unni
 */
public class CustomGeneInfoConverter extends BioFileConverter
{
    private static final String DATASET_TITLE = "Gene Info";
    private static final String DATA_SOURCE_NAME = "Gene Info";
    private String orgRefId;

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public CustomGeneInfoConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    private String geneSource = null;

    public void setGeneSource(String geneSource) {
        this.geneSource = geneSource;
    }

    public String getGeneSource() {
        return this.geneSource;
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        File currentFile = getCurrentFile();
        String taxonId = currentFile.getName().split("_")[0];
        try {
            // sanity check to ensure that the taxon ID is fetched
            long taxonIdInteger = Long.parseLong(taxonId);
        } catch (NumberFormatException e) {
            System.out.println("Unsupported file name encountered: " + currentFile.getName());
            System.out.println("The input file name should be of the format <TAXON_ID>_custom_gene_info.tab where <TAXON_ID> must be numerical (Ex: 9606_custom_gene_info.tab)");
            System.exit(1);
        }
        orgRefId = getOrganism(taxonId);

        Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            String primaryIdentifier = line[0].trim();
            //String secondaryIdentifier = line[1].trim();
            String geneSource = line[1].trim();
            String symbol = line[2].trim();
            String name = line[3].trim();
            Item geneItem = createItem("Gene");
            geneItem.setAttribute("primaryIdentifier", primaryIdentifier);
	    geneItem.setAttribute("source", geneSource);
            //if (!secondaryIdentifier.equals("-")) geneItem.setAttribute("secondaryIdentifier", secondaryIdentifier);
            if (!symbol.equals("-")) geneItem.setAttribute("symbol", symbol);
            if (!name.equals("-")) geneItem.setAttribute("name", name);
            geneItem.setReference("organism", orgRefId);
            try {
                store(geneItem);
            } catch(Exception e) {
                System.out.println("Error while storing item: " + geneItem + "\n" + "Exception stack trace:" + e);
            }
        }
    }
}
