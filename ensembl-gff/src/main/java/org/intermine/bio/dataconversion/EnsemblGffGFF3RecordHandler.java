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
 * A converter/retriever for the EnsemblGff dataset via GFF files.
 */

public class EnsemblGffGFF3RecordHandler extends GFF3RecordHandler
{
    /**
     * Create a new EnsemblGffGFF3RecordHandler for the given data model.
     * @param model the model for which items will be created
     */

    Map<String,String> aliasToRefId = new HashMap<String,String>();
    Map<String,String> geneToRefId = new HashMap<String,String>();
    Map<String,String> xRefToRefId = new HashMap<String,String>();

    public EnsemblGffGFF3RecordHandler (Model model) {
        super(model);
        refsAndCollections.put("Transcript", "gene");
        refsAndCollections.put("MRNA", "gene");
        refsAndCollections.put("TRNA", "gene");
        refsAndCollections.put("MiRNA", "gene");
        refsAndCollections.put("RRNA", "gene");
        refsAndCollections.put("SnRNA", "gene");
        refsAndCollections.put("SnoRNA", "gene");
        refsAndCollections.put("LincRNA", "gene");
        refsAndCollections.put("LncRNA", "gene");
        refsAndCollections.put("TRJGene", "gene");
        refsAndCollections.put("Ribozyme", "gene");
        refsAndCollections.put("ScaRNA", "gene");
        refsAndCollections.put("SRNA", "gene");
        refsAndCollections.put("TRCGene", "gene");
        refsAndCollections.put("TRVgene", "gene");
        refsAndCollections.put("IGCgene", "gene");
        refsAndCollections.put("IGVGene", "gene");
        refsAndCollections.put("Exon", "transcripts");
        refsAndCollections.put("CDS", "transcript");
        refsAndCollections.put("UTR", "transcripts");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(GFF3Record record) {

        Item feature = getFeature();
        String clsName = feature.getClassName();
        feature.removeAttribute("secondaryIdentifier");
        feature.removeAttribute("symbol");
        feature.setAttribute("source", record.getSource());
        if( clsName.equals("Gene") ) {
            if(record.getAttributes().get("ID") != null){
                String id = record.getAttributes().get("ID").iterator().next();
                feature.setAttribute("primaryIdentifier", id);
            }
            if(record.getAttributes().get("symbol_ensembl") != null){
                String symbol = record.getAttributes().get("symbol_ensembl").iterator().next();
                feature.setAttribute("symbol", symbol);
            }

            if (record.getAttributes().get("gene_biotype") != null) {
                String biotype = parseGeneBiotype(record.getAttributes().get("gene_biotype").iterator().next());
                feature.setAttribute("biotype", biotype);
            }

            if (record.getAttributes().get("description") != null) {
                String description = record.getAttributes().get("description").iterator().next();
                feature.setAttribute("description", URLDecoder.decode(URLDecoder.decode(description)));
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
        else if( clsName.equals("MRNA") || clsName.equals("SRNA") || clsName.equals("ScaRNA") ||clsName.equals("TRCGene") || clsName.equals("TRVGene") || clsName.equals("IGCGene") || clsName.equals("IGVGene") || clsName.equals("Ribozyme") || clsName.equals("Transcript") || clsName.equals("TRNA") || clsName.equals("MiRNA") || clsName.equals("RRNA") || clsName.equals("SnRNA") || clsName.equals("SnoRNA") || clsName.equals("LincRNA") || clsName.equals("LncRNA") || clsName.equals("TRJGene")) {
            if(record.getAttributes().get("ID") != null){
                String id = record.getAttributes().get("ID").iterator().next();
                feature.setAttribute("primaryIdentifier", id);
            }
            if(record.getAttributes().get("symbol_ensembl") != null){
                String symbol = record.getAttributes().get("symbol_ensembl").iterator().next();
                feature.setAttribute("symbol", symbol);
            }
            if(record.getAttributes().get("protein_id") != null){
                String prot = record.getAttributes().get("protein_id").iterator().next();
                feature.setAttribute("proteinIdentifier", prot);
            }
        }
        else if( clsName.equals("Exon") || clsName.equals("CDS") || clsName.equals("UTR") ) {
            // Do nothing
        }
        else {
            System.out.println("Unaccounted class type encountered: " + clsName);
            System.exit(1);
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
        if (biotype.equals("Mt_tRNA") || biotype.equals("TR_C_gene") || biotype.equals("TR_V_gene") ||biotype.equals("IG_C_gene") || biotype.equals("IG_V_gene")|| biotype.equals("Mt_rRNA") || biotype.equals("RNase_MRP_RNA") || biotype.equals("SRP_RNA") || biotype.equals("misc_RNA") || biotype.equals("C_region") || biotype.equals("V_segment") || biotype.equals("telomerase_RNA")) {
            returnType = biotype.replace("_", " ");
        }
        else if (biotype.equals("miRNA") || biotype.equals("sRNA") || biotype.equals("scaRNA") || biotype.equals("ribozyme") || biotype.equals("tRNA") || biotype.equals("rRNA") || biotype.equals("snRNA") || biotype.equals("snoRNA") || biotype.equals("lncRNA") || biotype.equals("lincRNA") || biotype.equals("lncRNA") || biotype.equals("TR_J_gene")) {
            returnType = biotype;
        }
        else if (biotype.equals("protein_coding") || biotype.equals("processed_pseudogene")) {
            String[] splitList = biotype.split("_");
            returnType = StringUtils.capitalize(splitList[0]) + " " + StringUtils.capitalize(splitList[1]);
        }
        else if (biotype.equals("pseudogene") || biotype.equals("other")) {
            returnType = StringUtils.capitalize(biotype);
        }
        else {
            System.out.println("Unexpected gene_biotype: " + biotype);
            System.exit(1);
        }
        return returnType;
    }
}
