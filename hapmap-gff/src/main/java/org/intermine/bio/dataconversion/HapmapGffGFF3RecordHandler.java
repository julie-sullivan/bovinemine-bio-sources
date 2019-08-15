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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.Iterator;

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
 * A converter/retriever for Hapmap GFF3
 */

public class HapmapGffGFF3RecordHandler extends GFF3RecordHandler
{
 Map<String,String> aliasToRefId = new HashMap<String,String>();

    /**
     * Create a new HapmapGFF3RecordHandler for the given data model.
     * @param model the model for which items will be created
     */
    public HapmapGffGFF3RecordHandler (Model model) {
        super(model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(GFF3Record record) {
        Item feature = getFeature();
        String clsName = feature.getClassName();
        feature.removeAttribute("symbol");

        if (clsName.equals("SNP")) {
            if (record.getAttributes().get("ID") != null) {
                String secondaryIdentifier = record.getAttributes().get("ID").iterator().next();
                feature.setAttribute("secondaryIdentifier", secondaryIdentifier);
            }
            if (record.getAttributes().get("Reference_seq") != null) {
                String refAllele = record.getAttributes().get("Reference_seq").iterator().next();
                feature.setAttribute("referenceAllele", refAllele);
            }
            if (record.getAttributes().get("Variant_seq") != null) {
                String altAllele = record.getAttributes().get("Variant_seq").iterator().next();
                feature.setAttribute("alternateAllele", altAllele);
            }

            List<String> dbxrefs = record.getDbxrefs();
            if (dbxrefs != null) {
                Iterator<String> dbxrefsIter = dbxrefs.iterator();
                while (dbxrefsIter.hasNext()) {
                    String dbxref = dbxrefsIter.next();
                    List<String> refList = new ArrayList<String>( Arrays.asList(StringUtil.split(dbxref, ",")));
                    for (String ref : refList) {
                        ref = ref.trim();
                        int colonIndex = ref.indexOf(":");
                        if (colonIndex == -1) {
                            throw new RuntimeException("Error in Dbxref attribute " + ref);
                        }
                        if (ref.startsWith("dbSNP")) {
                            feature.setAttribute("primaryIdentifier", ref.replace("dbSNP:", ""));
                        }
                    }
                }
            }

          if (record.getAliases() != null) {
                List<String> aliases = record.getAliases();
                Iterator<String> aliasesIterator = aliases.iterator();
                while (aliasesIterator.hasNext()) {
                    setAliasName(aliasesIterator.next());
                }
            }
        


       }

    }


       
     public void setAliasName(String alias) {
        // TODO: Should the relationship between Gene and AliasName be a Collection or a Reference?
        Item feature = getFeature();
        List<String> splitVal = new ArrayList<String>(Arrays.asList(StringUtil.split(alias, ":")));
        if (splitVal.size() != 2) {
            System.out.println("size: " + splitVal.size());
            System.out.println("Ambiguous aliasName: " + splitVal);
            System.out.println("Expected aliasName format is '<ALIAS_ID>:<ALIAS_SOURCE>'");
            System.out.println("Note: ALIAS_ID must be associated with its source");
            System.exit(1);
        }
        String aliasPrimaryIdentifier = splitVal.get(0);
        String aliasSource = splitVal.get(1);
            Item aliasItem = converter.createItem("AliasName");
            aliasItem.setAttribute("identifier", aliasPrimaryIdentifier);
            aliasItem.setAttribute("source", aliasSource);
            aliasItem.setReference("organism", getOrganism());
            String aliasRefId = aliasItem.getIdentifier();
            feature.addToCollection("alias", aliasRefId);
            aliasItem.addToCollection("gene", feature.getIdentifier());
            aliasToRefId.put(aliasPrimaryIdentifier, aliasRefId);
            addItem(aliasItem);
        
    }





}




