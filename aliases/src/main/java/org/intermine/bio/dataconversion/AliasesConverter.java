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
 * A Converter that loads aliases for Gene entity.
 * @author
 */
public class AliasesConverter extends BioFileConverter
{
    //
    private static final String DATASET_TITLE = null;
    private static final String DATA_SOURCE_NAME = null;
    Map<String, Item> geneItemMap = new HashMap<String, Item>();
    Map<String, Item> aliasItemMap = new HashMap<String, Item>();
    Map<String, Item> dataSourceItemMap = new HashMap<String, Item>();
    Map<String, Item> dataSetItemMap = new HashMap<String, Item>();
    String taxonId;
    String organismReferenceId;
    String currentAliasSource;
    Item dataSetItem;
    Item dataSourceItem;

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public AliasesConverter(ItemWriter writer, Model model) {
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
            String[] geneInfo = line[0].trim().split(":");
            String genePrimaryIdentifier = geneInfo[0];
            String geneSource = geneInfo[1];

            Item geneItem;
            if (geneItemMap.containsKey(genePrimaryIdentifier)) {
                geneItem = geneItemMap.get(genePrimaryIdentifier);
            }
            else {
                geneItem = createItem("Gene");
                geneItem.setAttribute("primaryIdentifier", genePrimaryIdentifier);
                geneItem.setAttribute("source", geneSource);
                geneItem.setReference("organism", organismReferenceId);
                geneItemMap.put(genePrimaryIdentifier, geneItem);
            }

            String[] aliasEntryList = line[1].trim().split(",");
            for (String aliasEntry : aliasEntryList) {
                String[] aliasInfo = aliasEntry.split(":");
                String str = aliasInfo[0];
                String info = str.replaceAll("@", ":");
                String aliasIdentifier = info;
                String aliasSource = aliasInfo[1];

                if (!aliasSource.equals(currentAliasSource)) {
                    dataSourceItem = getDataSourceItem(aliasSource);
                    dataSetItem = getDataSetItem(aliasSource, dataSourceItem);
                    currentAliasSource = aliasSource;
                }

                Item aliasItem;
                if (aliasItemMap.containsKey(aliasEntry)) {
                    aliasItem = aliasItemMap.get(aliasEntry);
                }
                else {
                    aliasItem = createItem("AliasName");
                    aliasItem.setAttribute("identifier", aliasIdentifier);
                    aliasItem.setAttribute("source", aliasSource);
                    aliasItem.setReference("organism", organismReferenceId);
                    aliasItem.addToCollection("dataSets", dataSetItem);
                }

                aliasItem.addToCollection("features", geneItem.getIdentifier());
                geneItem.addToCollection("aliases", aliasItem.getIdentifier());
                aliasItemMap.put(aliasEntry, aliasItem);
                geneItemMap.put(genePrimaryIdentifier, geneItem);
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
     * Get or create an item for a given dataSourceName
     * @param dataSourceName
     * @return
     */
    protected Item getDataSourceItem(String dataSourceName) {
        Item dataSourceItem;
        if (dataSourceName == null) {
            dataSourceItem = null;
        }
        else {
            if (dataSourceItemMap.containsKey(dataSourceName)) {
                dataSourceItem = dataSourceItemMap.get(dataSourceName);
            }
            else {
                dataSourceItem = createItem("DataSource");
                dataSourceItem.setAttribute("name", dataSourceName);
                dataSourceItemMap.put(dataSourceName, dataSourceItem);
            }
        }
        return dataSourceItem;
    }

    /**
     * Get or create an item for a given dataSetName
     * @param dataSetName
     * @param dataSourceItem
     * @return
     */
    protected Item getDataSetItem(String dataSetName, Item dataSourceItem) {
        Item dataSetItem;
        if (dataSetName == null) {
            dataSetItem = null;
        }
        else {
            if (dataSetItemMap.containsKey(dataSetName)) {
                dataSetItem = dataSetItemMap.get(dataSetName);
            }
            else {
                dataSetItem = createItem("DataSet");
                dataSetItem.setAttribute("name", dataSetName);
                dataSetItem.setReference("dataSource", dataSourceItem.getIdentifier());
                dataSetItemMap.put(dataSetName, dataSetItem);
            }
        }
        return dataSetItem;
    }

    /**
     *
     * {@inheritDoc}
     */
    @Override
    public void close() throws Exception {
        storeAllItems(dataSourceItemMap);
        storeAllItems(dataSetItemMap);
        storeAllItems(aliasItemMap);
        storeAllItems(geneItemMap);
    }
}
