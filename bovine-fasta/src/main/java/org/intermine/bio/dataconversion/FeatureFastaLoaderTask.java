package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2014 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.nbio.core.exceptions.ParserException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.DNASequenceCreator;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.io.PlainFastaHeaderParser;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.intermine.model.bio.Organism;
import org.intermine.objectstore.ObjectStoreException;

/**
 * Code for loading fasta for uniprot proteins.
 * @author julie
 */
public class FeatureFastaLoaderTask extends FastaLoaderTask
{
    private Map<String, Organism> organisms = new HashMap<String, Organism>();

    /**
     * {@inheritDoc}
     */
    @Override
    protected Organism getOrganism(Sequence bioJavaSequence) throws ObjectStoreException {
        String header = ((ProteinSequence) bioJavaSequence).getOriginalHeader();

        final String regexp = "OS\\=\\w+\\s\\w+";
        Pattern p = Pattern.compile(regexp);
        Matcher m = p.matcher(header);
        if (m.find()) {
            header = m.group();
            String[] bits = header.split("=");
            if (bits.length != 2) {
                return null;
            }

            String taxonId = getTaxonId(bits[1]);
            System.out.println(bits[1]);
            //System.exit(0);
            if (taxonId == null) {
                return null;
            }
            Organism org = organisms.get(taxonId);
            if (org == null) {
                System.out.println(taxonId);
                org = getDirectDataLoader().createObject(Organism.class);
                org.setTaxonId(taxonId);
                getDirectDataLoader().store(org);
                organisms.put(taxonId, org);
            }
            return org;
        }
        return null;
    }
}
