import biokt.CodonTables;
class CodonTableDump {
  public static void main(String[] args) {
    for (char a : "ACGT".toCharArray())
      for (char b : "ACGT".toCharArray())
        for (char c : "ACGT".toCharArray()) {
          String codon = new String(new char[] {a,b,c});
          System.out.println(codon + "\t" + CodonTables.INSTANCE.getSTANDARD().translate(codon));
        }
  }
}