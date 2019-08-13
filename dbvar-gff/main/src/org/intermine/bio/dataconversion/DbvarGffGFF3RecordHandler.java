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

import java.util.HashMap;
import java.util.HashSet;
import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.metadata.Model;
import org.intermine.metadata.StringUtil;
import org.intermine.xml.full.Item;
import java.util.Map;
import java.util.Map.Entry;

/**
 * A converter/retriever for the DbvarGff dataset via GFF files.
 */

public class DbvarGffGFF3RecordHandler extends GFF3RecordHandler
{

    /**
     * Create a new DbvarGffGFF3RecordHandler for the given data model.
     * @param model the model for which items will be created
     */
    public DbvarGffGFF3RecordHandler (Model model) {
        super(model);
        refsAndCollections.put("CopyNumberGain", "parent");
        refsAndCollections.put("CopyNumberLoss", "parent");
        refsAndCollections.put("Deletion", "parent");
        refsAndCollections.put("Inversion", "parent");
        refsAndCollections.put("TandemDuplicationVariant", "parent");
    }

    public HashMap<String, Item> topLevelItems =  new HashMap<String, Item>();
    public HashMap<String, Item> subItems =  new HashMap<String, Item>();

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(GFF3Record record) {
        // Even though the input file is a GVF, the format specifications are similar such that
        // GFF3Record can parse the data
        Item feature = getFeature();
        String clsName = feature.getClassName();
        feature.setAttribute("source", record.getSource());
        feature.removeAttribute("symbol");

        if (clsName.equals("CopyNumberVariation")) {
            String name = record.getAttributes().get("Name").iterator().next();
            // ID, by default, will be set as primaryIdentifier
            feature.setAttribute("name", name);
            if (record.getAttributes().get("validated") != null) {
                ArrayList<String> validated = new ArrayList<String> (record.getAttributes().get("validated"));
                if (validated.size() > 1) {
                    feature.setAttribute("validationStatus", StringUtil.join(validated, ","));
                }
                else {
                    feature.setAttribute("validationStatus", record.getAttributes().get("validated").iterator().next());
                }
            }
            if (record.getAttributes().get("sample_name") != null) {
                feature.setAttribute("sampleName", record.getAttributes().get("sample_name").get(0));
            }
            if (record.getAttributes().get("gender") != null) {
                feature.setAttribute("gender", record.getAttributes().get("gender").get(0));
            }
            if (record.getAttributes().get("remapScore") != null) {
                feature.setAttribute("remapScore", record.getAttributes().get("remapScore").iterator().next());
            }
            if (record.getAttributes().get("phenotype") != null) {
                feature.setAttribute("phenotype", record.getAttributes().get("phenotype").iterator().next());
            }
            if (record.getAttributes().get("isComplexSubstitution") != null) {
                feature.setAttribute("isComplexSubstitution", record.getAttributes().get("isComplexSubstitution").iterator().next());
            }
            if (record.getAttributes().get("Start_range") != null) {
                String innerStart = record.getAttributes().get("Start_range").get(0);
                String outerStart = record.getAttributes().get("Start_range").get(1);
                if (!innerStart.equals(".")) {
                    feature.setAttribute("innerStartCoordinate", innerStart);
                }
                if (!outerStart.equals(".")) {
                    feature.setAttribute("outerStartCoordinate", outerStart);
                }
            }
            if (record.getAttributes().get("End_range") != null) {
                String innerEnd = record.getAttributes().get("End_range").get(0);
                String outerEnd = record.getAttributes().get("End_range").get(1);
                if (!innerEnd.equals(".")) {
                    feature.setAttribute("innerStartCoordinate", innerEnd);
                }
                if (!outerEnd.equals(".")) {
                    feature.setAttribute("outerStartCoordinate", outerEnd);
                }
            }
        }
        else if (clsName.equals("CopyNumberGain") || clsName.equals("CopyNumberLoss")) {
            feature.setAttribute("name", record.getAttributes().get("Name").iterator().next());
            if (record.getAttributes().get("validated") != null) {
                ArrayList<String> validated = new ArrayList<String> (record.getAttributes().get("validated"));
                if (validated.size() > 1) {
                    feature.setAttribute("validationStatus", StringUtil.join(validated, ","));
                }
                else {
                    feature.setAttribute("validationStatus", record.getAttributes().get("validated").iterator().next());
                }
            }
            if (record.getAttributes().get("sample_name") != null) {
                feature.setAttribute("sampleName", record.getAttributes().get("sample_name").get(0));
            }
            if (record.getAttributes().get("gender") != null) {
                feature.setAttribute("gender", record.getAttributes().get("gender").iterator().next());
            }
            if (record.getAttributes().get("remapScore") != null) {
                feature.setAttribute("remapScore", record.getAttributes().get("remapScore").iterator().next());
            }
            if (record.getAttributes().get("phenotype") != null) {
                feature.setAttribute("phenotype", record.getAttributes().get("phenotype").iterator().next());
            }
            if (record.getAttributes().get("copy_number") != null) {
                feature.setAttribute("copyNumber", record.getAttributes().get("copy_number").iterator().next());
            }
            if (record.getAttributes().get("Start_range") != null) {
                String innerStart = record.getAttributes().get("Start_range").get(0);
                String outerStart = record.getAttributes().get("Start_range").get(1);
                if (!innerStart.equals(".")) {
                    feature.setAttribute("innerStartCoordinate", innerStart);
                }
                if (!outerStart.equals(".")) {
                    feature.setAttribute("outerStartCoordinate", outerStart);
                }
            }
            if (record.getAttributes().get("End_range") != null) {
                String innerEnd = record.getAttributes().get("End_range").get(0);
                String outerEnd = record.getAttributes().get("End_range").get(1);
                if (!innerEnd.equals(".")) {
                    feature.setAttribute("innerStopCoordinate", innerEnd);
                }
                if (!outerEnd.equals(".")) {
                    feature.setAttribute("outerStopCoordinate", outerEnd);
                }
            }
        }
        else if (clsName.equals("Deletion")) {
            feature.setAttribute("name", record.getAttributes().get("Name").iterator().next());
            if (record.getAttributes().get("validated") != null) {
                ArrayList<String> validated = new ArrayList<String> (record.getAttributes().get("validated"));
                if (validated.size() > 1) {
                    feature.setAttribute("validationStatus", StringUtil.join(validated, ","));
                }
                else {
                    feature.setAttribute("validationStatus", record.getAttributes().get("validated").iterator().next());
                }
            }
            if (record.getAttributes().get("sample_name") != null) {
                feature.setAttribute("sampleName", record.getAttributes().get("sample_name").get(0));
            }
            if (record.getAttributes().get("gender") != null) {
                feature.setAttribute("gender", record.getAttributes().get("gender").iterator().next());
            }
            if (record.getAttributes().get("remapScore") != null) {
                feature.setAttribute("remapScore", record.getAttributes().get("remapScore").iterator().next());
            }
            if (record.getAttributes().get("phenotype") != null) {
                feature.setAttribute("phenotype", record.getAttributes().get("phenotype").iterator().next());
            }
            if (record.getAttributes().get("Start_range") != null) {
                String innerStart = record.getAttributes().get("Start_range").get(0);
                String outerStart = record.getAttributes().get("Start_range").get(1);
                if (!innerStart.equals(".")) {
                    feature.setAttribute("innerStartCoordinate", innerStart);
                }
                if (!outerStart.equals(".")) {
                    feature.setAttribute("outerStartCoordinate", outerStart);
                }
            }
            if (record.getAttributes().get("End_range") != null) {
                String innerEnd = record.getAttributes().get("End_range").get(0);
                String outerEnd = record.getAttributes().get("End_range").get(1);
                if (!innerEnd.equals(".")) {
                    feature.setAttribute("innerStopCoordinate", innerEnd);
                }
                if (!outerEnd.equals(".")) {
                    feature.setAttribute("outerStopCoordinate", outerEnd);
                }
            }
        }
        else if (clsName.equals("Inversion")) {
            feature.setAttribute("name", record.getAttributes().get("Name").iterator().next());
            if (record.getAttributes().get("validated") != null) {
                ArrayList<String> validated = new ArrayList<String> (record.getAttributes().get("validated"));
                if (validated.size() > 1) {
                    feature.setAttribute("validationStatus", StringUtil.join(validated, ","));
                }
                else {
                    feature.setAttribute("validationStatus", record.getAttributes().get("validated").iterator().next());
                }
            }
            if (record.getAttributes().get("sample_name") != null) {
                feature.setAttribute("sampleName", record.getAttributes().get("sample_name").get(0));
            }
            if (record.getAttributes().get("gender") != null) {
                feature.setAttribute("gender", record.getAttributes().get("gender").iterator().next());
            }
            if (record.getAttributes().get("remapScore") != null) {
                feature.setAttribute("remapScore", record.getAttributes().get("remapScore").iterator().next());
            }
            if (record.getAttributes().get("phenotype") != null) {
                feature.setAttribute("phenotype", record.getAttributes().get("phenotype").iterator().next());
            }
            if (record.getAttributes().get("copy_number") != null) {
                feature.setAttribute("copyNumber", record.getAttributes().get("copy_number").iterator().next());
            }
            if (record.getAttributes().get("Start_range") != null) {
                String innerStart = record.getAttributes().get("Start_range").get(0);
                String outerStart = record.getAttributes().get("Start_range").get(1);
                if (!innerStart.equals(".")) {
                    feature.setAttribute("innerStartCoordinate", innerStart);
                }
                if (!outerStart.equals(".")) {
                    feature.setAttribute("outerStartCoordinate", outerStart);
                }
            }
            if (record.getAttributes().get("End_range") != null) {
                String innerEnd = record.getAttributes().get("End_range").get(0);
                String outerEnd = record.getAttributes().get("End_range").get(1);
                if (!innerEnd.equals(".")) {
                    feature.setAttribute("innerStopCoordinate", innerEnd);
                }
                if (!outerEnd.equals(".")) {
                    feature.setAttribute("outerStopCoordinate", outerEnd);
                }
            }
        }
        else if (clsName.equals("TandemDuplication")) {
            feature.setAttribute("name", record.getAttributes().get("Name").iterator().next());
            if (record.getAttributes().get("validated") != null) {
                ArrayList<String> validated = new ArrayList<String> (record.getAttributes().get("validated"));
                if (validated.size() > 1) {
                    feature.setAttribute("validationStatus", StringUtil.join(validated, ","));
                }
                else {
                    feature.setAttribute("validationStatus", record.getAttributes().get("validated").iterator().next());
                }
            }
            if (record.getAttributes().get("sample_name") != null) {
                feature.setAttribute("sampleName", record.getAttributes().get("sample_name").get(0));
            }
            if (record.getAttributes().get("gender") != null) {
                feature.setAttribute("gender", record.getAttributes().get("gender").iterator().next());
            }
            if (record.getAttributes().get("remapScore") != null) {
                feature.setAttribute("remapScore", record.getAttributes().get("remapScore").iterator().next());
            }
            if (record.getAttributes().get("phenotype") != null) {
                feature.setAttribute("phenotype", record.getAttributes().get("phenotype").iterator().next());
            }
            if (record.getAttributes().get("copy_number") != null) {
                feature.setAttribute("copyNumber", record.getAttributes().get("copy_number").iterator().next());
            }
            if (record.getAttributes().get("Start_range") != null) {
                String innerStart = record.getAttributes().get("Start_range").get(0);
                String outerStart = record.getAttributes().get("Start_range").get(1);
                if (!innerStart.equals(".")) {
                    feature.setAttribute("innerStartCoordinate", innerStart);
                }
                if (!outerStart.equals(".")) {
                    feature.setAttribute("outerStartCoordinate", outerStart);
                }
            }
            if (record.getAttributes().get("End_range") != null) {
                String innerEnd = record.getAttributes().get("End_range").get(0);
                String outerEnd = record.getAttributes().get("End_range").get(1);
                if (!innerEnd.equals(".")) {
                    feature.setAttribute("innerStopCoordinate", innerEnd);
                }
                if (!outerEnd.equals(".")) {
                    feature.setAttribute("outerStopCoordinate", outerEnd);
                }
            }
        }
        else if (clsName.equals("TandemDuplicationVariant")) {
            feature.setAttribute("name", record.getAttributes().get("Name").iterator().next());
            feature.setAttribute("source", record.getSource());
            if (record.getAttributes().get("validated") != null) {
                ArrayList<String> validated = new ArrayList<String> (record.getAttributes().get("validated"));
                if (validated.size() > 1) {
                    feature.setAttribute("validationStatus", StringUtil.join(validated, ","));
                }
                else {
                    feature.setAttribute("validationStatus", record.getAttributes().get("validated").iterator().next());
                }
            }
            if (record.getAttributes().get("sample_name") != null) {
                feature.setAttribute("sampleName", record.getAttributes().get("sample_name").get(0));
            }
            if (record.getAttributes().get("gender") != null) {
                feature.setAttribute("gender", record.getAttributes().get("gender").iterator().next());
            }
            if (record.getAttributes().get("remapScore") != null) {
                feature.setAttribute("remapScore", record.getAttributes().get("remapScore").iterator().next());
            }
            if (record.getAttributes().get("phenotype") != null) {
                feature.setAttribute("phenotype", record.getAttributes().get("phenotype").iterator().next());
            }
            if (record.getAttributes().get("copy_number") != null) {
                feature.setAttribute("copyNumber", record.getAttributes().get("copy_number").iterator().next());
            }
            if (record.getAttributes().get("Start_range") != null) {
                String innerStart = record.getAttributes().get("Start_range").get(0);
                String outerStart = record.getAttributes().get("Start_range").get(1);
                if (!innerStart.equals(".")) {
                    feature.setAttribute("innerStartCoordinate", innerStart);
                }
                if (!outerStart.equals(".")) {
                    feature.setAttribute("outerStartCoordinate", outerStart);
                }
            }
            if (record.getAttributes().get("End_range") != null) {
                String innerEnd = record.getAttributes().get("End_range").get(0);
                String outerEnd = record.getAttributes().get("End_range").get(1);
                if (!innerEnd.equals(".")) {
                    feature.setAttribute("innerStopCoordinate", innerEnd);
                }
                if (!outerEnd.equals(".")) {
                    feature.setAttribute("outerStopCoordinate", outerEnd);
                }
            }
        }
        else {
            System.out.println("Unaccounted feature type: " + clsName);
            System.exit(1);
        }
    }

}
