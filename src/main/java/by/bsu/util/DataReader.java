package by.bsu.util;

import by.bsu.model.IlluminaSNVSample;
import by.bsu.model.PairEndRead;
import by.bsu.model.Sample;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Class to read input data from standard input files
 */
public class DataReader {

    public static String[] readList(Path filePath) throws IOException {
        return readList(filePath, false);
    }

    public static String[] readList(Path filePath, boolean reversed) throws IOException {
        List<String> raw = Files.readAllLines(filePath);
        List<String> result = new ArrayList<>();
        int i = 0;
        StringBuilder seq = new StringBuilder();
        for (int j = 0; j < raw.size(); j++) {
            String str = raw.get(j);
            if (str.startsWith(">") && seq.length() > 0 || str.length() == 0) {
                result.add((reversed ? seq.reverse() : seq).toString());
                i++;
                seq.setLength(0);
            } else if (!str.startsWith(">")) {
                seq.append(str);
                //workaround for for cases when file doesn't contain last empty string
                if (j + 1 == raw.size()) {
                    result.add((reversed ? seq.reverse() : seq).toString());
                }
            }
        }
        return result.stream().toArray(String[]::new);
    }

    public static String[] readList(String filePath) throws IOException {
        return readList(Paths.get(filePath));
    }

    public static List<Sample> readSampleList(String filePath, boolean onePerFile) throws IOException {
        return readSampleList(filePath, onePerFile, false);
    }

    public static List<Sample> readSampleList(String filePath, boolean onePerFile, boolean reversed) throws IOException {
        File dir = new File(filePath);
        List<Sample> result = new ArrayList<>();
        for (File file : dir.listFiles()) {
            if (onePerFile == file.isDirectory()) {
                continue;
            }
            result.add(onePerFile ? readSampleFromFile(file) : readSampleFromFolder(file, reversed));
        }
        return result;
    }

    public static Sample readSampleFromFolder(File file) {
        return readSampleFromFolder(file, false);
    }

    public static Sample readSampleFromFolder(File file, boolean reversed) {
        List<File> sortedFiles = Arrays.stream(file.listFiles())
                .sorted(Comparator.comparing(File::getName))
                .collect(Collectors.toList());
        List<String> seq = new ArrayList<>();
        sortedFiles.stream().filter(f -> !f.isHidden()).forEach(f -> {
            try {
                seq.addAll(Arrays.asList(readList(f.toPath(), reversed)));
            } catch (IOException e) {
                System.err.println("Error reading file for multyfile sample");
                e.printStackTrace();
            }
        });

        return new Sample(file.getName(), seq.stream().toArray(String[]::new));
    }

    public static Sample readSampleFromFile(File file) throws IOException {
        return new Sample(file.getName(), readList(file.toPath()));
    }

    public static Sample readOneLined(File file) throws IOException {
        List<String> raw = Files.readAllLines(file.toPath());
        List<String> tmp = new ArrayList<>();
        List<String> N = new ArrayList<>();
        List<Integer> offset = new ArrayList<>();
        List<Integer> readLength = new ArrayList<>();
        N.add("");
        for (int i = 1; i < 3000; i++) {
            N.add(N.get(i - 1) + 'N');
        }
        int i = 0;
        for (String s : raw) {
            String[] split = s.split("\t");
            offset.add(Integer.valueOf(split[3]) - 1);
            readLength.add(split[5].length());
            tmp.add(N.get(offset.get(i)) + split[5]);
            i++;
        }

        return new Sample(file.getName(),
                Utils.stringsForHamming(tmp.toArray(new String[0]), 'N'));
    }

    /**
     * Splits each column of the given sample into 4 columns for each minor. For example:
     * Given alphabet as "ACGT-"
     * for i-th column 'G' is major. then i-th column will split into 4 columns where first column will place '2' for each 'A' and '1' for any other allele,
     * for second column will place '2' for each 'C' and '1' otherwise, for third '1' for each 'T' and for forth '1' for each '-'
     *
     * @param sample
     * @param profile
     * @param alphabet
     * @return
     */
    public static Sample splitColumns(Sample sample, double[][] profile, String alphabet) {
        String[] sequences = new String[sample.sequences.length];
        int i = 0;
        for (String sequence : sample.sequences) {
            sequences[i++] = getSplittedRead(sequence, 0, profile, alphabet);
        }
        return new Sample(sample.name + "_splitted", sequences);
    }

    public static IlluminaSNVSample splitColumns(IlluminaSNVSample sample, double[][] profile, String alphabet) {
        int minorCount = alphabet.length() - 1;
        List<PairEndRead> splittedReads = new ArrayList<>();
        sample.reads.forEach(r -> {
            splittedReads.add(new PairEndRead(getSplittedRead(r.l, r.lOffset, profile, alphabet), getSplittedRead(r.r, r.rOffset, profile, alphabet),
                    r.lOffset * minorCount, r.rOffset * minorCount, r.name));
        });
        return new IlluminaSNVSample(sample.name + "_splitted", splittedReads, sample.referenceLength * minorCount);
    }

    public static IlluminaSNVSample getIlluminaPairedReads(File file) {
        Map<String, List<SAMRecord>> readsSet = new HashMap<>();
        System.out.println("Start read sam");
        SamReader open = SamReaderFactory.make().open(file);
        for (SAMRecord anOpen : open) {
            if (readsSet.size() % 1_000 == 0) {
                System.out.print("\r" + readsSet.size());
            }
            if (anOpen.getReadUnmappedFlag()) {
                continue;
            }
            if (!readsSet.containsKey(anOpen.getReadName())) {
                readsSet.put(anOpen.getReadName(), new ArrayList<>());
            }
            readsSet.get(anOpen.getReadName()).add(anOpen);
        }
        System.out.println(" DONE");
        SAM4WebLogo sam4WebLogo = new SAM4WebLogo(open);
        System.out.println("Start convert");
        List<PairEndRead> pairedReads = new ArrayList<>(readsSet.size());
        Pattern begin = Pattern.compile("^-*");
        Pattern end = Pattern.compile("-*$");
        readsSet.forEach((key, value) -> {
            if (value.size() == 1) {
                //TODO rewrite this shit
                String s = sam4WebLogo.printRead(value.get(0));
                s = end.matcher(begin.matcher(s).replaceAll("")).replaceAll("");
                pairedReads.add(new PairEndRead(s,
                        "",
                        value.get(0).getAlignmentStart() - 1,
                        -1, key)
                );
            } else {
                for (int i = 0; i < value.size(); i++) {
                    for (int j = 0; j < value.size(); j++) {
                        if (value.get(i).getMateAlignmentStart() == value.get(j).getAlignmentStart()
                                && value.get(i).getAlignmentStart() <= value.get(j).getAlignmentStart()) {
                            if (i == j) {
                                String s = sam4WebLogo.printRead(value.get(i));
                                s = end.matcher(begin.matcher(s).replaceAll("")).replaceAll("");
                                pairedReads.add(new PairEndRead(s,
                                        "",
                                        value.get(i).getAlignmentStart() - 1,
                                        -1, key)
                                );
                            } else {
                                String s = sam4WebLogo.printRead(value.get(i));
                                s = end.matcher(begin.matcher(s).replaceAll("")).replaceAll("");
                                String s2 = sam4WebLogo.printRead(value.get(j));
                                s2 = end.matcher(begin.matcher(s2).replaceAll("")).replaceAll("");
                                pairedReads.add(new PairEndRead(s,
                                        s2,
                                        value.get(i).getAlignmentStart() - 1,
                                        value.get(j).getAlignmentStart() - 1,
                                        key)
                                );
                            }
                        }
                    }
                }
                //TODO think about it later
//                for (int i = 0; i < l.getBaseQualities().length; i++) {
//                    if (l.getBaseQualities()[i] < 30){
//                        strL.replace(i, i+1, "N");
//                    }
//                }
            }
            if (pairedReads.size() % 1_000 == 0) {
                System.out.print("\r" + pairedReads.size());
            }
        });
        System.out.println(" DONE");
        return new IlluminaSNVSample(file.getName(), pairedReads,
                pairedReads.parallelStream().mapToInt(r -> Math.max(r.lOffset + r.l.length(), r.rOffset + r.r.length())).max().orElse(0));
    }

    private static String getSplittedRead(String read, int offset, double[][] profile, String alphabet) {
        StringBuilder str = new StringBuilder();
        for (int j = 0; j < read.length(); j++) {
            int major = Utils.getMajorAllele(profile, offset + j);
            int minor = 0;
            for (int k = 0; k < alphabet.length() - 1; k++, minor++) {
                if (minor == major) {
                    minor++;
                }
                int allele = alphabet.indexOf(read.charAt(j));
                str.append(allele == minor ? "2" : "1");
            }
        }
        return str.toString();
    }
}
