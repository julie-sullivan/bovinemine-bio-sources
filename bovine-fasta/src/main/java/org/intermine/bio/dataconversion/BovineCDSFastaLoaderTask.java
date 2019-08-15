package org.intermine.bio.dataconversion;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.intermine.metadata.Model;
import org.intermine.model.FastPathObject;
import org.intermine.model.InterMineObject;
import org.intermine.objectstore.ObjectStore;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.DynamicUtil;

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


/**
 * A fasta loader that understand the headers of CDS fasta files and can make the
 * appropriate extra objects and references.
 * @author Kim Rutherford
 * This script has been adapted with the help of AIPCDSFastaLoaderTask.java and FlyBaseCDSFastaLoaderTask.java
 */
public class BovineCDSFastaLoaderTask extends BovineFeatureFastaLoaderTask
{
   /**
     * {@inheritDoc}
     */
    @Override
    protected void extraProcessing (Sequence bioJavaSequence,
                                    org.intermine.model.bio.Sequence flymineSequence,
                                    BioEntity bioEntity, Organism organism, DataSet dataSet)
    {
//        String mrnaIdentifier = bioJavaSequence.getAccession().getID();
//        bioJavaSequence.
//        String regexp = "^.+\\s+(\\S+):([0-9]+-[0-9]+)\\s+(\\S+)\\s+(\\S+)\\s+hasEarlyStopCodon=(.+)$";
//        Pattern p = Pattern.compile(regexp);
//        Matcher m = p.matcher(header);
//        String hasESC = "";
//        if (m.matches()) {
//            hasESC = m.group(5);
//        }
//        if (hasESC != "") {
//            bioEntity.setFieldValue("hasEarlyStopCodon", hasESC.toLowerCase());
//        }
//        ObjectStore os = getIntegrationWriter().getObjectStore();
//        Model model = os.getModel();
//        if (model.hasClassDescriptor(model.getPackageName() + ".CDS")) {
//            Class<? extends FastPathObject> cdsCls =
//                model.getClassDescriptorByName("CDS").getType();
//            if (!DynamicUtil.isInstance(bioEntity, cdsCls)) {
//                throw new RuntimeException("the InterMineObject passed to "
//                        + "BovineCDSFastaLoaderTask.extraProcessing() is not a "
//                        + "CDS: " + bioEntity);
//            }
//            InterMineObject mrna = getMRNA(mrnaIdentifier, getGeneSource(), organism, model);
//            if (mrna != null) {
//                bioEntity.setFieldValue("transcript", mrna);
//            }
//
//            Location loc = getLocationFromHeader(header, (SequenceFeature) bioEntity,
//                    organism);
//            getDirectDataLoader().store(loc);
//        } else {
//            throw new RuntimeException("Trying to load CDS sequence but CDS does not exist in the"
//                    + " data model");
//        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected String getIdentifier(Sequence bioJavaSequence) {
//        Annotation annotation = bioJavaSequence.getAnnotation();
//        String mrnaIdentifier = bioJavaSequence.getName();
//        String header = (String) annotation.getProperty("description");
//        String last = header.substring(header.lastIndexOf(' ') + 1);
//        return mrnaIdentifier + "-CDS";

        // TODO this is wrong
        return bioJavaSequence.getAccession().getID();
    }
}
