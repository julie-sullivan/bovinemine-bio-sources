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

import java.lang.System;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.net.URLDecoder;
import org.apache.commons.lang.StringUtils;

import java.util.HashMap;
import java.util.HashSet;
import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.metadata.Model;
import org.intermine.metadata.StringUtil;
import org.intermine.xml.full.Item;
import java.util.Map;
import java.util.Map.Entry;

/**
 * A converter/retriever for the RefseqNoncodingGff dataset via GFF files.
 */

public class RefseqNoncodingGffGFF3RecordHandler extends GFF3RecordHandler
{

    /**
     * Create a new RefseqNoncodingGffGFF3RecordHandler for the given data model.
     * @param model the model for which items will be created
     */
    Map<String,String> aliasToRefId = new HashMap<String,String>();
    Map<String,String> xRefToRefId = new HashMap<String,String>();
    Map<String,String> geneToRefId = new HashMap<String,String>();

    public RefseqNoncodingGffGFF3RecordHandler (Model model) {
        super(model);
        refsAndCollections.put("Exon", "transcripts");
        refsAndCollections.put("Transcript", "gene");
        refsAndCollections.put("PrimaryTranscript", "gene");
        refsAndCollections.put("NcRNA", "gene");
        refsAndCollections.put("MiRNA", "gene");
        refsAndCollections.put("TRNA", "gene");
        refsAndCollections.put("RRNA", "gene");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(GFF3Record record) {
        Item feature = getFeature();
        String clsName = feature.getClassName();
        feature.setAttribute("source", record.getSource());

        if (clsName.equals("Gene")) {

            if (record.getAttributes().get("symbol_ncbi") != null) {
                String symbol = record.getAttributes().get("symbol_ncbi").iterator().next();
                feature.setAttribute("symbol", symbol);
            }
            else {
                feature.removeAttribute("symbol");
            }

            if (record.getAttributes().get("gene_biotype") != null) {
                String biotype = parseGeneBiotype(record.getAttributes().get("gene_biotype").iterator().next());
                feature.setAttribute("biotype", biotype);
            }

            if (record.getAttributes().get("description") != null) {
                String description = record.getAttributes().get("description").iterator().next();
                feature.setAttribute("description", URLDecoder.decode(URLDecoder.decode(description)));
            }

            if (record.getAttributes().get("duplicate_entity") != null) {
                String duplicates = record.getAttributes().get("duplicate_entity").iterator().next();
                List<String> entities = new ArrayList<String>( Arrays.asList(StringUtil.split(duplicates, "|")));
                for (String entity : entities) {
                    Item duplicateEntityItem = converter.createItem("DuplicateEntity");
                    String duplicateEntityItemRefId = duplicateEntityItem.getIdentifier();
                    List<String> entityAttributes = new ArrayList<String>( Arrays.asList(StringUtil.split(entity, ",")));
                    List<String> locationInformation = new ArrayList<String>( Arrays.asList(StringUtil.split(entityAttributes.get(0), ":")));
                    String chromosome = locationInformation.get(0);
                    List<String> positionInfo = new ArrayList<String>( Arrays.asList(StringUtil.split(locationInformation.get(1), "..")));
                    String start = positionInfo.get(0);
                    String end = positionInfo.get(1);
                    int strand = locationInformation.get(2).equals("+") ? 1 : -1;
                    duplicateEntityItem.setAttribute("chromosome", chromosome);
                    duplicateEntityItem.setAttribute("start", start);
                    duplicateEntityItem.setAttribute("end", end);
                    duplicateEntityItem.setAttribute("strand", Integer.toString(strand));

                    if (entityAttributes.size() > 1) {
                        String geneIdentifier = entityAttributes.get(1).replace("NCBI_Gene:", "");
                        duplicateEntityItem.setAttribute("geneIdentifier", geneIdentifier);
                    }
                    if (entityAttributes.size() > 2) {
                        ArrayList<String> product = new ArrayList<String>(Arrays.asList(StringUtil.split(entityAttributes.get(2), ":")));

                        if (product.size() > 1) {
                            String transcriptId = product.get(1);
                            duplicateEntityItem.setAttribute("transcriptIdentifier", transcriptId);
                        }
                        if (product.size() > 2) {
                            String proteinIdentifier = product.get(2);
                            duplicateEntityItem.setAttribute("proteinIdentifier", proteinIdentifier);
                        }
                    }
                    try {
                        converter.store(duplicateEntityItem);
                    }
                    catch (Exception e) {
                        System.out.println("Exception while storing duplicateEntityItem:" + duplicateEntityItem + "\n" + e);
                    }
                    feature.addToCollection("duplicateEntities", duplicateEntityItemRefId);
                }
            }

            // Accessing Dbxrefs
            List<String> dbxrefs = record.getDbxrefs();
            if (dbxrefs != null) {
                Iterator<String> dbxrefsIter = dbxrefs.iterator();
                while (dbxrefsIter.hasNext()) {
                    String dbxref = dbxrefsIter.next();
                    List<String> refList = new ArrayList<String>( Arrays.asList(StringUtil.split(dbxref, ",")));
                    for( String ref : refList) {
                        ref = ref.trim();
                        int colonIndex = ref.indexOf(":");
                        if( colonIndex == -1 ) {
                            throw new RuntimeException("Error in Dbxref attribute " + ref );
                        }
                        if( ref.startsWith("NCBI_Gene") ) {
                            feature.setAttribute( "primaryIdentifier", ref.replace("NCBI_Gene:", "") );
                        }
                        if( ref.startsWith("BGD:") ) {
                            // treating as xRef
                            // TODO: Verify this assumption
                            //String xRefValue = ref.replace("BGD:", "") +  ":btau_OGSv2";
                            //setCrossReference(xRefValue);
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

            if (record.getAttributes().get("xRef") != null) {
                List<String> xRefList = record.getAttributes().get("xRef");
                Iterator<String> xRefIterator = xRefList.iterator();
                while (xRefIterator.hasNext()) {
                    setCrossReference(xRefIterator.next());
                }
            }
        }
        else if (clsName.equals("Transcript")) {

            if (record.getAttributes().get("symbol_ncbi") != null) {
                String symbol = record.getAttributes().get("symbol_ncbi").iterator().next();
                feature.setAttribute("symbol", symbol);
            }
            else {
                feature.removeAttribute("symbol");
            }
            if (record.getAttributes().get("description") != null) {
                String description = record.getAttributes().get("description").iterator().next();
                feature.setAttribute("description", URLDecoder.decode(URLDecoder.decode(description)));
            }

//            if (record.getAttributes().get("feature_type") != null) {
//                String ft = record.getAttributes().get("feature_type").iterator().next();
//                feature.setAttribute("status", ft);
//            }

            // Accessing Dbxrefs
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
                            throw new RuntimeException("Error in Dbxref attribute " + ref );
                        }

                        if (ref.startsWith("RefSeq_NA:")) {
                            feature.setAttribute("primaryIdentifier", ref.replace("RefSeq_NA:", ""));
                        }
                    }
                }
            }
        }
        else if (clsName.equals("PrimaryTranscript")) {

            if (record.getAttributes().get("symbol_ncbi") != null) {
                String symbol = record.getAttributes().get("symbol_ncbi").iterator().next();
                feature.setAttribute("symbol", symbol);
            }
            else {
                feature.removeAttribute("symbol");
            }
            if (record.getAttributes().get("description") != null) {
                String description = record.getAttributes().get("description").iterator().next();
                feature.setAttribute("description", URLDecoder.decode(URLDecoder.decode(description)));
            }

//            if (record.getAttributes().get("feature_type") != null) {
//                String ft = record.getAttributes().get("feature_type").iterator().next();
//                feature.setAttribute("status", ft);
//            }

            if (record.getAttributes().get("mature_form") != null) {
                // mature_form for miRNA
                String matureFormString = record.getAttributes().get("mature_form").iterator().next();
                List<String> entities = new ArrayList<String>(Arrays.asList(StringUtil.split(matureFormString, "|")));
                for (String entity : entities) {
                    Item matureTranscriptItem = converter.createItem("MatureTranscript");
                    String matureTranscriptItemRefId = matureTranscriptItem.getIdentifier();
                    List<String> entityAttributes = new ArrayList<String>(Arrays.asList(StringUtil.split(entity, ",")));
                    List<String> locationInformation = new ArrayList<String>(Arrays.asList(StringUtil.split(entityAttributes.get(0), ":")));
                    String chromosome = locationInformation.get(0);
                    List<String> positionInfo = new ArrayList<String>(Arrays.asList(StringUtil.split(locationInformation.get(1), "..")));
                    String start = positionInfo.get(0);
                    String end = positionInfo.get(1);
                    int strand = locationInformation.get(2).equals("+") ? 1 : -1;
                    matureTranscriptItem.setAttribute("chromosome", chromosome);
                    matureTranscriptItem.setAttribute("start", start);
                    matureTranscriptItem.setAttribute("end", end);
                    matureTranscriptItem.setAttribute("strand", Integer.toString(strand));
                    String transcriptId = entityAttributes.get(1).split(":")[1];
                    String mirbaseId = entityAttributes.get(2).split(":")[1];
                    String description = entityAttributes.get(3).split(":")[1];
                    matureTranscriptItem.setAttribute("transcriptIdentifier", transcriptId);
                    matureTranscriptItem.setAttribute("mirbaseIdentifier", mirbaseId);
                    matureTranscriptItem.setAttribute("description", URLDecoder.decode(URLDecoder.decode(description)));

                    try {
                        converter.store(matureTranscriptItem);
                    }
                    catch (Exception e) {
                        System.out.println("Exception while storing duplicateEntityItem:" + matureTranscriptItem + "\n" + e);
                    }
                    feature.addToCollection("matureTranscripts", matureTranscriptItemRefId);
                }
            }

            // Accessing Dbxrefs
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
                            throw new RuntimeException("Error in Dbxref attribute " + ref );
                        }
                        if (ref.startsWith("miRBase:")) {
                            feature.setAttribute("mirbaseIdentifier", ref.replace("miRBase:", ""));
                        }
                        if (ref.startsWith("RefSeq_NA:")) {
                            feature.setAttribute("primaryIdentifier", ref.replace("RefSeq_NA:", ""));
                        }
                    }
                }
            }
        }
        else if (clsName.equals("MiRNA")) {

            if (record.getAttributes().get("symbol_ncbi") != null) {
                String symbol = record.getAttributes().get("symbol_ncbi").iterator().next();
                feature.setAttribute("symbol", symbol);
            }
            else {
                feature.removeAttribute("symbol");
            }
            if (record.getAttributes().get("description") != null) {
                String description = record.getAttributes().get("description").iterator().next();
                feature.setAttribute("description", URLDecoder.decode(URLDecoder.decode(description)));
            }

//            if (record.getAttributes().get("feature_type") != null) {
//                String ft = record.getAttributes().get("feature_type").iterator().next();
//                feature.setAttribute("status", ft);
//            }

            // Accessing Dbxrefs
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
                            throw new RuntimeException("Error in Dbxref attribute " + ref );
                        }
                        if (ref.startsWith("miRBase")) {
                            feature.setAttribute("primaryIdentifier", ref.replace("miRBase:", ""));
                        }
                    }
                }
            }
        }
        else if (clsName.equals("TRNA")) {

            if (record.getAttributes().get("symbol_ncbi") != null) {
                String symbol = record.getAttributes().get("symbol_ncbi").iterator().next();
                feature.setAttribute("symbol", symbol);
            }
            else {
                feature.removeAttribute("symbol");
            }
            if (record.getAttributes().get("description") != null) {
                String description = record.getAttributes().get("description").iterator().next();
                feature.setAttribute("description", URLDecoder.decode(URLDecoder.decode(description)));
            }

//            if (record.getAttributes().get("feature_type") != null) {
//                String ft = record.getAttributes().get("feature_type").iterator().next();
//                feature.setAttribute("status", ft);
//            }

            // Accessing Dbxrefs
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
                            throw new RuntimeException("Error in Dbxref attribute " + ref );
                        }
                        if (ref.startsWith("RefSeq_NA")) {
                            feature.setAttribute("primaryIdentifier", ref.replace("RefSeq_NA:", ""));
                        }
                    }
                }
            }
        }
        else if (clsName.equals("RRNA")) {

            if( record.getAttributes().get("symbol_ncbi") != null) {
                String symbol = record.getAttributes().get("symbol_ncbi").iterator().next();
                feature.setAttribute("symbol", symbol);
            }
            else {
                feature.removeAttribute("symbol");
            }
            if (record.getAttributes().get("description") != null) {
                String description = record.getAttributes().get("description").iterator().next();
                feature.setAttribute("description", URLDecoder.decode(URLDecoder.decode(description)));
            }

//            if (record.getAttributes().get("feature_type") != null) {
//                String ft = record.getAttributes().get("feature_type").iterator().next();
//                feature.setAttribute("status", ft);
//            }

            // Accessing Dbxrefs
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
                            throw new RuntimeException("Error in Dbxref attribute " + ref );
                        }
                        if (ref.startsWith("RefSeq_NA")) {
                            feature.setAttribute("primaryIdentifier", ref.replace("RefSeq_NA:", ""));
                        }
                    }
                }
            }
        }
        else if( clsName.equals("Exon") ) {
            // Do nothing
        }
        else {
            System.out.println("Unexpected feature type encountered: " + clsName);
            System.exit(0);
        }
    }

    /**
     * Method parses the alias string, creates an AliasName item and sets the necessary references and collections
     * @param alias
     */
    public void setAliasName(String alias) {
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
        if (aliasToRefId.containsKey(aliasPrimaryIdentifier)) {
            feature.addToCollection("aliases", aliasToRefId.get(aliasPrimaryIdentifier));
        } else {
            Item aliasItem = converter.createItem("AliasName");
            aliasItem.setAttribute("identifier", aliasPrimaryIdentifier);
            aliasItem.setAttribute("source", aliasSource);
            aliasItem.setReference("organism", getOrganism());
            String aliasRefId = aliasItem.getIdentifier();
            feature.addToCollection("aliases", aliasRefId);
            aliasItem.addToCollection("features", feature.getIdentifier());
            aliasToRefId.put(aliasPrimaryIdentifier, aliasRefId);
            addItem(aliasItem);
        }
    }

    /**
     * Method parses the xRef string, creates a xRef item, creates a Gene item and sets the necessary references and collections
     * @param xRef
     */
    public void setCrossReference(String xRef) {
        Item feature = getFeature();
        List<String> xRefPair = new ArrayList<String>(Arrays.asList(StringUtil.split(xRef, ":")));
        if (xRefPair.size() == 0) { return; }
        if (xRefPair.size() != 2) {
            System.out.println("Ambiguous xRef: " + xRefPair);
            System.out.println("Expected xRef format is '<XREF_ID>:<XREF_SOURCE>'");
            System.out.println("Note: XREF_SOURCE should match column 2 of the alternate GFF3 (if any)");
            System.exit(1);
        }
        String identifier = xRefPair.get(0);
        String xRefSource = xRefPair.get(1);
        if (xRefToRefId.containsKey(identifier)) {
            feature.addToCollection("dbCrossReferences", xRefToRefId.get(identifier));
            if (! geneToRefId.containsKey(identifier)) {
                System.out.println("xRef exists but its corresponding gene instance does not exist");
                System.exit(1);
            }
        } else {
            Item xRefItem = converter.createItem("xRef");
            xRefItem.setAttribute("refereeSource", xRefSource);
            xRefItem.setReference("organism", getOrganism());
            String xRefRefId = xRefItem.getIdentifier();
            feature.addToCollection("dbCrossReferences", xRefRefId);
            xRefToRefId.put(identifier, xRefRefId);
            if (!geneToRefId.containsKey(identifier)) {
                // storing the Gene instance of xRef
                Item geneItem = converter.createItem("Gene");
                geneItem.setAttribute("primaryIdentifier", identifier);
                geneItem.setAttribute("source", xRefSource);
                geneItem.setReference("organism", getOrganism());
                geneToRefId.put(identifier, geneItem.getIdentifier());
                xRefItem.setReference("referrer", feature.getIdentifier());
                xRefItem.setReference("referee", geneItem.getIdentifier());
                addItem(geneItem);
            }
            addItem(xRefItem);
        }
    }

    /**
     * Parse biotype from gene_biotype attribute for aesthetics
     * @param biotype
     * @return
     */
    public String parseGeneBiotype(String biotype) {
        String returnType = "";
        if (biotype.equals("guide_RNA") || biotype.equals("antisense_RNA") || biotype.equals("Mt_tRNA") || biotype.equals("Mt_rRNA") || biotype.equals("RNase_MRP_RNA") || biotype.equals("SRP_RNA") || biotype.equals("misc_RNA") || biotype.equals("C_region") || biotype.equals("V_segment") || biotype.equals("telomerase_RNA")) {
            returnType = biotype.replace("_", " ");
        }
        else if (biotype.equals("miRNA") || biotype.equals("tRNA") || biotype.equals("rRNA") || biotype.equals("snRNA") || biotype.equals("snoRNA") || biotype.equals("lncRNA")) {
            returnType = biotype;
        }
        else if (biotype.equals("protein_coding") || biotype.equals("processed_pseudogene")) {
            String[] splitList = biotype.split("_");
            returnType = StringUtils.capitalize(splitList[0]) + " " + StringUtils.capitalize(splitList[1]);
        }
        else if (biotype.equals("pseudogene") || biotype.equals("transcribed_pseudogene") || biotype.equals("other")) {
            returnType = StringUtils.capitalize(biotype);
        }
        else {
            System.out.println("Unexpected gene_biotype: " + biotype);
            System.exit(1);
        }
        return returnType;
    }
}
