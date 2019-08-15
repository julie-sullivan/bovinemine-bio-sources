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
import org.intermine.xml.full.Item;
//import org.apache.commons.lang3.StringUtils;
import java.util.List;
/**
 * A converter/retriever for the OgsGff dataset via GFF files.
 */

public class OgsGffGFF3RecordHandler extends GFF3RecordHandler
{

    /**
     * Create a new OgsGffGFF3RecordHandler for the given data model.
     * @param model the model for which items will be created
     */
    public OgsGffGFF3RecordHandler (Model model) {
        super(model);
        refsAndCollections.put("MRNA", "gene");
        refsAndCollections.put("Exon", "transcripts");
        refsAndCollections.put("CDS", "transcript");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(GFF3Record record) {
        Item feature = getFeature();
        String clsName = feature.getClassName();
        feature.setAttribute("source", record.getSource());
        feature.removeAttribute("symbol");

        if( clsName.equals("Gene") ) {
            if(record.getAttributes().get("ID") != null){
                String id = record.getAttributes().get("ID").iterator().next();
                feature.setAttribute("primaryIdentifier", id);
            }
            if(record.getAttributes().get("Name") != null){
                String name = record.getAttributes().get("Name").iterator().next();
                feature.setAttribute("name", name);
            }
            if(record.getAttributes().get("description") != null){
                //String description = record.getAttributes().get("description").iterator().next();
                List<String> description = record.getAttributes().get("description");
	
                //feature.setAttribute("description", description);
                feature.setAttribute("description", String.join(",", description));
            }
            if(record.getAttributes().get("gene_biotype") != null){
                String biotype = record.getAttributes().get("gene_biotype").iterator().next();
                feature.setAttribute("biotype", biotype);
            }
        }
        else if( clsName.equals("MRNA") || clsName.equals("Polypeptide") ) {
            if(record.getAttributes().get("ID") != null){
                String id = record.getAttributes().get("ID").iterator().next();
                feature.setAttribute("primaryIdentifier", id);
            }
            if(record.getAttributes().get("description") != null){
                //String description = record.getAttributes().get("description").iterator().next();
                List<String> description = record.getAttributes().get("description");
                //feature.setAttribute("description", description);
                feature.setAttribute("description", String.join(",", description));
            }
            if(record.getAttributes().get("transcript_biotype") != null){
                String biotype = record.getAttributes().get("transcript_biotype").iterator().next();
                feature.setAttribute("biotype", biotype);
            }
            if(record.getAttributes().get("source") != null){
                String source = record.getAttributes().get("source").iterator().next();
            }
            if(record.getAttributes().get("Name") != null){
                String name = record.getAttributes().get("Name").iterator().next();
                feature.setAttribute("name", name);
            }
        }
        else if( clsName.equals("StartCodon") || clsName.equals("StopCodon") ) {
            // Do nothing
        }
    }
}
