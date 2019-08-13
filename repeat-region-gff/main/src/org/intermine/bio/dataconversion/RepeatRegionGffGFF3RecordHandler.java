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
 * A converter/retriever for the RepeatRegionGff dataset via GFF files.
 */

public class RepeatRegionGffGFF3RecordHandler extends GFF3RecordHandler
{

    /**
     * Create a new RepeatRegionGffGFF3RecordHandler for the given data model.
     * @param model the model for which items will be created
     */
    public RepeatRegionGffGFF3RecordHandler (Model model) {
        super(model);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(GFF3Record record) {
        Item feature = getFeature();
        String clsName = feature.getClassName();
        feature.setAttribute("source", record.getSource());

        if (clsName.equals("RepeatRegion")) {
            String target = record.getAttributes().get("target").iterator().next();
            String name = target.split(" ")[0];
            feature.setAttribute("name", name);
            feature.setAttribute("primaryIdentifier", name.replace("Motif:", ""));
            feature.removeAttribute("symbol");
        }
    }

    @Override
    public void setLocation(Item location) {
        location.setAttribute("doNotComputeOverlaps", "Y");
        items.put("_location", location);
    }
}
