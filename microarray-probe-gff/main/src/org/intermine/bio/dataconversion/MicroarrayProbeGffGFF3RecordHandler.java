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

/**
 * A converter/retriever for the MicroarrayProbeGff dataset via GFF files.
 */

public class MicroarrayProbeGffGFF3RecordHandler extends GFF3RecordHandler
{

    /**
     * Create a new MicroarrayProbeGffGFF3RecordHandler for the given data model.
     * @param model the model for which items will be created
     */
    public MicroarrayProbeGffGFF3RecordHandler (Model model) {
        super(model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(GFF3Record record) {
        Item feature = getFeature();
        String clsName = feature.getClassName();
        if (clsName.equals("Probe")) {
            if (record.getAttributes().get("Name") != null) {
                String name = record.getAttributes().get("Name").iterator().next();
                feature.setAttribute("name", name);
                feature.setAttribute("primaryIdentifier", name);
            }
            feature.setAttribute("source", record.getSource());
            feature.removeAttribute("symbol");
        }
    }
}
