package edu.caltech.lncrna.arraytools.programs;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.caltech.lncrna.bio.alignment.Aligned;
import edu.caltech.lncrna.bio.alignment.Alignment;
import edu.caltech.lncrna.bio.annotation.Annotated;
import edu.caltech.lncrna.bio.annotation.Annotation;
import edu.caltech.lncrna.bio.annotation.BedFileRecord;
import edu.caltech.lncrna.bio.datastructures.GenomeTree;
import edu.caltech.lncrna.bio.io.BamParser;
import edu.caltech.lncrna.bio.io.BedParser;

/**
This program was originally written to assist Joanna with designing RAP probes
that target transposable elements. In typical probe design, we discard probes
that align to multiple locations in the genome. Transposons, though, have many
copies scattered throughout the genome, and we need to know the nature (type of
element, whether in a gene, etc) of every location that each probe aligns to.
 */
public class TransposonProbeAnalyzer {

    private final Path repeatsPath;
    private final Path genesPath;
    private final Path probesPath;
    
    private final GenomeTree<BedFileRecord> genes;
    private final GenomeTree<Annotated> introns;
    private final GenomeTree<BedFileRecord> repeats;
    private final Map<String, Probe> probes;
    
    private static final String VERSION = "1.0.0";
    private static final Logger LOGGER = Logger.getLogger("TransposonProbeAnalyzer");
    
    private static final String HELP_TEXT = "java -jar "
            + "TransposonProbeAnalyzer.jar --genes genes.bed --probes "
            + "probes.bam --repeats repeats.bed > output.txt";
    
    public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        CommandLine cmd = parseArgs(args);
        TransposonProbeAnalyzer program = new TransposonProbeAnalyzer(cmd);
        program.loadRepeats();
        program.loadGenes();
        program.loadProbes();
        program.print();
        LOGGER.info("Program complete");
        LOGGER.info((System.currentTimeMillis() - startTime) + " milliseconds elapsed.");
    }
    
    public TransposonProbeAnalyzer(CommandLine cmd) {
        
        genes = new GenomeTree<>();
        introns = new GenomeTree<>();
        repeats = new GenomeTree<>();
        probes = new HashMap<>();
        
        repeatsPath = Paths.get(cmd.getOptionValue("repeats"));
        genesPath = Paths.get(cmd.getOptionValue("genes"));
        probesPath = Paths.get(cmd.getOptionValue("probes"));
        
        if (cmd.hasOption("debug")) {
            LOGGER.setLevel(Level.FINEST);
            LOGGER.info("Running in debug mode.");
        }

    }
    
    private static CommandLine parseArgs(String[] args) {
        
        Option versionOption = Option.builder("v")
                .longOpt("version")
                .desc("show version information")
                .hasArg(false)
                .required(false)
                .build();
        
        Option helpOption = Option.builder("h")
                .longOpt("help")
                .desc("show usage information")
                .hasArg(false)
                .required(false)
                .build();
        
        Option debugOption = Option.builder()
                .longOpt("debug")
                .desc("debug mode")
                .hasArg(false)
                .required(false)
                .build();
        
        Option repeatsOption = Option.builder()
                .longOpt("repeats")
                .desc("the BED file of repeat regions")
                .hasArg(true)
                .required(true)
                .build();
        
        Option genesOption = Option.builder()
                .longOpt("genes")
                .desc("the BED file of gene regions")
                .hasArg(true)
                .required(true)
                .build();
        
        Option probesOption = Option.builder()
                .longOpt("probes")
                .desc("the BAM file of probe alignments")
                .hasArg(true)
                .required(true)
                .build();
        
        Options helpOptions = new Options().addOption(helpOption);
        Options versionOptions = new Options().addOption(versionOption);
        
        Options mainOptions = new Options()
                .addOption(repeatsOption)
                .addOption(genesOption)
                .addOption(probesOption)
                .addOption(debugOption);
        
        Options allOptions = new Options();
        helpOptions.getOptions().forEach(allOptions::addOption);
        versionOptions.getOptions().forEach(allOptions::addOption);
        mainOptions.getOptions().forEach(allOptions::addOption);
        HelpFormatter formatter = new HelpFormatter();
                
        CommandLineParser parser = new DefaultParser();
        CommandLine rtrn = null;
        
        try {
            CommandLine cmds = parser.parse(helpOptions, args, true);
            if (cmds.getOptions().length == 1) {
                formatter.printHelp(HELP_TEXT, allOptions);
                System.exit(0);
            }
            cmds = parser.parse(versionOptions, args, true);
            if (cmds.getOptions().length == 1) {
                System.out.println(VERSION);
                System.exit(0);
            }
            rtrn = parser.parse(mainOptions, args);
        } catch (ParseException e) {
            formatter.printHelp(HELP_TEXT, allOptions);
            System.exit(1);
        }
        
        if (rtrn == null) {
            LOGGER.severe("An unknown error occurred while parsing command line arguments");
            System.exit(1);
        }
        
        return rtrn;
    }
    
    public void loadRepeats() {
        LOGGER.info("Loading repeats.");
        try (BedParser bp = new BedParser(repeatsPath)) {
            bp.stream().forEach(repeats::add);
        }
        LOGGER.info("Loaded " + repeats.size() + " repeat annotations");
    }
    
    public void loadGenes() {
        LOGGER.info("Loading genes.");
        try (BedParser bp = new BedParser(genesPath)) {
            bp.stream().forEach(genes::add);
        }
        LOGGER.info("Loaded " + genes.size() + " gene annotations.");
        
        LOGGER.info("Loading introns.");
        for (Annotated a : genes) {
            a.getIntronStream().forEach(introns::add);
        }
        LOGGER.info("Loaded " + introns.size() + " intron annotations.");
    }
    
    public void loadProbes() {
        LOGGER.info("Loading probes.");
        try (BamParser<? extends Aligned<? extends Alignment>> bp =
                BamParser.newInstance(probesPath)) {
            bp.getAlignmentStream()
                .peek(x -> LOGGER.finest("Reading probe " + x.toFormattedBedString(4)))
                .forEach(x -> addRead(x));
        }
        LOGGER.info("Loaded " + probes.size() + " probes.");
    }
    
    public void addRead(Alignment a) {
        LOGGER.log(Level.FINEST, "Adding probe " + a.getName());
        String name = a.getName();
        Probe probe = probes.getOrDefault(name, new Probe());
        probe.addPosition(a);
        probes.put(name, probe);
    }
    
    public void print() {
        System.out.println("NAME\tREPEATS\tGENES_NO_REPEATS\t" + 
                "GENES_EXONS_WITH_REPEATS\tGENES_INTRONS_WITH_REPEATS");
        for (Probe probe : probes.values()) {
            System.out.println(probe.toString());
        }
    }
    
    /**
     * A class to represent a RAP probe.
     */
    public class Probe {
        
        private String name;
        
        /**
         * A list of positions that this probe aligns to.
         */
        private final List<Position> positions;
        
        public Probe() {
            positions = new ArrayList<>();
        }
        
        public void addPosition(Alignment a) {
            name = a.getName();
            positions.add(new Position(a));
        }
        
        /**
         * A string representation of this probe suitable for printing into the
         * output text file.
         */
        @Override
        public String toString() {
            
            Map<String, Integer> repeatNames = new HashMap<>();
            Map<String, Integer> geneNoRepeatsNames = new HashMap<>();
            Map<String, Integer> geneExonOverlapsRepeatNames = new HashMap<>();
            Map<String, Integer> geneIntronOverlapsRepeatNames = new HashMap<>();

            for (Position position : positions) {
                
                for (String geneName : position.geneNoRepeatsNames) {
                    int val = geneNoRepeatsNames.getOrDefault(geneName, 0);
                    geneNoRepeatsNames.put(geneName, val + 1);
                }
                
                for (String geneName : position.geneExonOverlapsRepeatNames) {
                    int val = geneExonOverlapsRepeatNames.getOrDefault(geneName, 0);
                    geneExonOverlapsRepeatNames.put(geneName, val + 1);
                }
                
                for (String geneName : position.geneIntronOverlapsRepeatNames) {
                    int val = geneIntronOverlapsRepeatNames.getOrDefault(geneName, 0);
                    geneIntronOverlapsRepeatNames.put(geneName, val + 1);
                }
                
                for (String repeatName : position.repeatNames) {
                    int val = repeatNames.getOrDefault(repeatName, 0);
                    repeatNames.put(repeatName, val + 1);
                }
            }
            
            StringBuilder sb = new StringBuilder();
            sb.append(name + "\t");

            if (repeatNames.isEmpty()) {
                sb.append(".");
            } else {
                for (Entry<String, Integer> entry : repeatNames.entrySet()) {
                    sb.append(entry.getValue() + ":" + entry.getKey() + ";");
                }
            }
            
            sb.append("\t");
            
            if (geneNoRepeatsNames.isEmpty()) {
                sb.append(".");
            } else {
                for (Entry<String, Integer> entry : geneNoRepeatsNames.entrySet()) {
                    sb.append(entry.getValue() + ":" + entry.getKey() + ";");
                }
            }
            
            sb.append("\t");

            if (geneExonOverlapsRepeatNames.isEmpty()) {
                sb.append(".");
            } else {
                for (Entry<String, Integer> entry : geneExonOverlapsRepeatNames.entrySet()) {
                    sb.append(entry.getValue() + ":" + entry.getKey() + ";");
                }
            }
            
            sb.append("\t");

            if (geneIntronOverlapsRepeatNames.isEmpty()) {
                sb.append(".");
            } else {
                for (Entry<String, Integer> entry : geneIntronOverlapsRepeatNames.entrySet()) {
                    sb.append(entry.getValue() + ":" + entry.getKey() + ";");
                }
            }
            
            return sb.toString();
        }
        
        @Override
        public boolean equals(Object other) {
            if (this == other) {
                return true;
            }
            
            if (!(other instanceof Probe)) {
                return false;
            }
            
            Probe o = (Probe) other;
            return positions.equals(o.positions);
        }
        
        @Override
        public int hashCode() {
            return positions.hashCode();
        }
    }
    
    /**
     * An extension of the {@link Annotation} class that keeps track of the
     * names of any overlapping elements or genes. Used to represent where a
     * probe aligns.
     */
    public class Position extends Annotation {
        
        private final String name;
        private final List<String> geneNoRepeatsNames;
        private final List<String> geneExonOverlapsRepeatNames;
        private final List<String> geneIntronOverlapsRepeatNames;
        private final List<String> repeatNames;
        
        public Position(Alignment a) {
            super(a);
            name = a.getName();
            geneNoRepeatsNames = new ArrayList<>();
            geneExonOverlapsRepeatNames = new ArrayList<>();
            geneIntronOverlapsRepeatNames = new ArrayList<>();
            repeatNames = new ArrayList<>();
            repeats.bodyOverlappers(a.getBody())
                .forEachRemaining(x -> repeatNames.add(x.getName()));
            
            genes.bodyOverlappers(a.getBody())
                 .forEachRemaining(x -> {
                     if (!repeats.overlaps(x.getBody())) {
                         geneNoRepeatsNames.add(x.getName());
                     } else if (repeats.overlaps(x)) {
                         geneExonOverlapsRepeatNames.add(x.getName());
                     } else {
                         geneIntronOverlapsRepeatNames.add(x.getName());
                     }
                 });
        }
        
        @Override
        public boolean equals(Object other) {
            
            if (this == other) {
                return true;
            }
            
            if (!(other instanceof Position)) {
                return false;
            }
            
            Position o = (Position) other;
            return super.equals(o) &&
                    name.equals(o.name) &&
                    geneNoRepeatsNames.equals(o.geneNoRepeatsNames) &&
                    geneExonOverlapsRepeatNames.equals(o.geneExonOverlapsRepeatNames) &&
                    geneIntronOverlapsRepeatNames.equals(o.geneIntronOverlapsRepeatNames) &&
                    repeatNames.equals(o.repeatNames);
        }
        
        @Override
        public int hashCode() {
            int hashCode = super.hashCode();
            hashCode = 37 * hashCode + name.hashCode();
            hashCode = 37 * hashCode + geneNoRepeatsNames.hashCode();
            hashCode = 37 * hashCode + geneExonOverlapsRepeatNames.hashCode();
            hashCode = 37 * hashCode + geneIntronOverlapsRepeatNames.hashCode();
            hashCode = 37 * hashCode + repeatNames.hashCode();
            return hashCode;
        }
    }
}
