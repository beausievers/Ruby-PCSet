#
# => PCSet Class for Ruby
# => Beau Sievers
# => Hanover, Fall 2008.
#
#
# TODO: Make this a module to avoid namespace collisions.
#       Lilypond and MusicXML output
#

include Math

def choose(n, k)
  return [[]] if n.nil? || n.empty? && k == 0
  return [] if n.nil? || n.empty? && k > 0
  return [[]] if n.size > 0 && k == 0
  c2 = n.clone
  c2.pop
  new_element = n.clone.pop
  choose(c2, k) + append_all(choose(c2, k-1), new_element)
end

def append_all(lists, element)
  lists.map { |l| l << element }
end

def arrayToBinary(array)
  array.inject(0) {|sum, n| sum + 2**n}
end

# the following method is horrifically inelegant 
# but avoids monkey-patching.
# TODO: do this right, incl. good error checking
def pearsons(x, y)
  if !x.is_a?(Array) || !y.is_a?(Array) then raise StandardError, "x and y must be arrays", caller end
  if x.size != y.size then raise StandardError, "x and y must be same size", caller end
  sumX = x.inject(0) {|sum, n| sum + n}
  sumY = y.inject(0) {|sum, n| sum + n}
  sumSquareX = x.inject(0) {|sum, n| sum + n * n}
  sumSquareY = y.inject(0) {|sum, n| sum + n * n}
  xy = []
  x.zip(y) {|a, b| xy.push(a * b)}
  sumXY = xy.inject(0) {|sum, n| sum + n}
  num = sumXY - ((sumX * sumY)/x.size)
  den = Math.sqrt((sumSquareX - ((sumX*sumX)/x.size)) * (sumSquareY - ((sumY*sumY)/x.size)))
  (num/den)
end

class PCSet
  include Comparable
  attr_reader :pitches, :base, :input
  
  def initialize(pcarray, base=12)
    if pcarray.instance_of?(Array) && pcarray.all?{|pc| pc.instance_of?(Fixnum)}
      @base, @input = base, pcarray
      @pitches = pcarray.map{ |x| x % @base }.uniq
    else
      raise ArgumentError, "Improperly formatted PC array", caller
    end
  end
  
  def PCSet.newFromString(pcstring, base=12)
    if base > 36 then raise StandardError, "Use PCSet.new to create pcsets with a base larger than 36", caller end
    pcarray = []
    pcstring.downcase.split(//).each do |c|
      if c <= 'z' and c >= '0' then pcarray.push(c.to_i(36)) end
    end
    PCSet.new pcarray, base
  end
  
  def <=>(pcs)
    @pitches <=> pcs.pitches
  end
  
  def [](index)
    @pitches[index]
  end

  # Intersection
  def &(other)
    PCSet.new @pitches & other.pitches
  end
  
  # Union
  def |(other)
    PCSet.new @pitches | other.pitches
  end
  
  def inspect
    @pitches.inspect
  end
  
  def length  
    @pitches.length
  end
  
  def invert(axis=0)
    PCSet.new @pitches.map {|x| (axis-x) % @base}
  end
  
  def invert!(axis=0)
    @pitches.map! {|x| (axis-x) % @base}
  end
    
  def transpose(interval)
    PCSet.new @pitches.map {|x| (x + interval) % @base}
  end
  
  def transpose!(interval)
    @pitches.map! {|x| (x + interval) % @base}
  end
  
  def multiply(m = 5)
    PCSet.new @pitches.map {|x| (x * m) % @base}
  end
  
  def multiply!(m = 5)
    @pitches.map! {|x| (x * m) % @base}
  end
  
  def zero
    transpose(-1 * @pitches[0])
  end
  
  def zero!
    transpose!(-1 * @pitches[0])
  end
  
  def transpositions
    (0..(@base-1)).to_a.map{|x| @pitches.map {|y| (y + x) % @base}}.sort.map {|x| PCSet.new x}
  end
  
  def transpositionsAndInversions(axis=0)
    transpositions + invert(axis).transpositions
  end
  
  #
  # Normal form after Straus. Morris and AthenaCL do this differently.
  #
  def normalForm
    tempar = @pitches.sort
    arar = []                 # [[1,4,7,8,10],[4,7,8,10,1], etc.] get each cyclic variation
    tempar.each {arar.push PCSet.new(tempar.unshift(tempar.pop))}
    mostLeftCompact(arar)
  end
  
  def normalForm!
    @pitches = normalForm.pitches
  end
  
  def isNormalForm?
    self.pitches == self.normalForm.pitches
  end
  
  def setClass
    transpositionsAndInversions.map{|pcs| pcs.normalForm}.sort
  end
  
  def prime
    mostLeftCompact([normalForm.zero, invert.normalForm.zero])
  end
  
  def prime!
    self.pitches = self.prime.pitches
  end
  
  def isPrime?
    self.pitches == self.prime.pitches
  end
  
  def complement
    newPitches = []
    @base.times do |p|
      if !@pitches.include? p then
        newPitches.push p
      end
    end
    PCSet.new newPitches
  end
  
  def fullIntervalVector
    pairs = choose(@pitches, 2)                       # choose every pc pair
    intervals = pairs.map {|x| (x[1] - x[0]) % @base} # calculate every interval
    iVector = Array.new(@base-1).fill(0)
    intervals.each {|x| iVector[x-1] += 1}            # count the intervals
    iVector
  end
  
  def intervalVector
    iVector = fullIntervalVector
    (0..((@base-1)/2)-1).each {|x| iVector[x] += iVector.pop}
    iVector
  end
  
  #
  # Morris's invariance vector
  #
  def invarianceVector(m = 5)
    t = transpositions.map!{|pcs| self & pcs}
    ti = invert.transpositions.map!{|pcs| self & pcs}    
    tm = multiply(m).transpositions.map!{|pcs| self & pcs}   
    tmi = invert.multiply(m).transpositions.map!{|pcs| self & pcs}   
    tc = complement.transpositions.map!{|pcs| self & pcs}   
    tic = complement.invert.transpositions.map!{|pcs| self & pcs}
    tmc = complement.multiply(m).transpositions.map!{|pcs| self & pcs}
    tmic = complement.invert.multiply(m).transpositions.map!{|pcs| self & pcs}
    
    [t, ti, tm, tmi, tc, tic, tmc, tmic].map{|x| x.reject{|pcs| pcs.pitches != @pitches}.length}
  end

  # Huron's aggregate dyadic consonance measure. 
  # Huron. Interval-Class Content in Equally Tempered Pitch-Class Sets: 
  # Common Scales Exhibit Optimum Tonal Consonance. 
  # Music Perception (1994) vol. 11 (3) pp. 289-305
  def huron
    if @base != 12 then raise StandardError, "PCSet.huron only makes sense for mod 12 pcsets", caller end
      
    #              m2/M7   M2/m7  m3/M6  M3/m6  P4/P5   A4/d5
    huronTable = [-1.428, -0.582, 0.594, 0.386, 1.240, -0.453]
    intervalConsonance = []
    intervalVector.zip(huronTable) {|x, y| intervalConsonance.push(x * y) }
    aggregateDyadicConsonance = intervalConsonance.inject {|sum, n| sum + n}
    [aggregateDyadicConsonance, pearsons(intervalVector, huronTable)]
  end

  #
  # Balzano's vector of relations. Citation for all Balzano methods:
  # 
  # Balzano. "The Pitch Set as a Level of Description for Studying Musical
  # Pitch Perception" in Music, Mind, and Brain ed. Clynes. Plenum Press. 1982.
  #
  def vectorOfRelations
    (0..length-1).to_a.map do |i|
      (0..length-1).to_a.map do |j|
        (@pitches[(i + j) % length] - @pitches[i]) % @base
      end
    end
  end

  #
  # Checks if the set satisfies Balzano's uniqueness.
  #
  def isUnique?
    vectorOfRelations.uniq.size == vectorOfRelations.size
  end
  
  #
  # Checks if the set satisfies Balzano's scalestep-semitone coherence.
  # For all s[i] and s[i1]:
  # j < k => v[i][j] < v[i1][k]
  # Where j and k are scalestep-counting indices.
  # And unless v[i][j] == 6 (a tritone), in which case the strict inequality is relaxed.
  #
  def isCoherent?
    v = vectorOfRelations
    truthArray = []
    allPairIndices = choose((0..length-1).to_a, 2)
    allPairIndices.each do |i, i1|
      allPairIndices.each do |j, k|
        if v[i][j] == 6
          truthArray.push(v[i][j] <= v[i1][k])
        else
          truthArray.push(v[i][j] < v[i1][k])    
        end
        if v[i1][j] == 6
          truthArray.push(v[i1][j] <= v[i][k])
        else
          truthArray.push(v[i1][j] < v[i][k])
        end
      end
    end
    !truthArray.include?(false)
  end
  
  #
  # Strict Balzano coherence, no inequality relaxation for tritones.
  #
  def isStrictlyCoherent?
    v = vectorOfRelations
    truthArray = []
    allPairIndices = choose((0..length-1).to_a, 2)
    allPairIndices.each do |i, i1|
      allPairIndices.each do |j, k|
        truthArray.push(v[i][j] < v[i1][k])
        truthArray.push(v[i1][j] < v[i][k])
      end
    end
    !truthArray.include?(false)
  end
  
  def notes(middleC=0)
    noteArray = ['C','C#','D','D#','E','F','F#','G','G#','A','A#','B']
    if @base != 12 then raise StandardError, "PCSet.notes only makes sense for mod 12 pcsets", caller end
    outString = String.new
    transpose(-middleC).pitches.each do |p|
      outString += noteArray[p] + ", "
    end
    outString.chop.chop
  end

  def info
    print "modulo: #{@base}\n"
    print "raw input: #{@input.inspect}\n"
    print "pitch set: #{@pitches.inspect}\n"
    print "notes: #{notes}\n"
    print "normal: #{normalForm.inspect}\n"
    print "prime: #{prime.inspect}\n"
    print "interval vector: #{intervalVector.inspect}\n"
    print "invariance vector: #{invarianceVector.inspect}\n"
    print "huron ADC: #{huron[0]}  pearsons: #{huron[1]}\n"
    print "balzano coherence: "
    if isStrictlyCoherent?
      print "strictly coherent\n"
    elsif isCoherent?
      print "coherent\n"
    else
      print "false\n"
    end
  end
  
  
#  def lilypond
#  
#  end
#  
#  def musicXML
#    
#  end

###############################################################################
  private
  
  #
  # Convert every pitch array to a binary representation, e.g.:
  # [0,2,4,8,10] -> 010100010101
  #           2^n:  BA9876543210
  # The smallest binary number is the most left-compact.
  #
  def mostLeftCompact(pcSetArray)
    if !pcSetArray.all? {|pcs| pcs.length == pcSetArray[0].length}
      raise ArgumentError, "PCSet.mostLeftCompact: All PCSets must be of same cardinality", caller
    end
    zeroedPitchArrays = pcSetArray.map {|pcs| pcs.zero.pitches}
    binaries = zeroedPitchArrays.map {|array| arrayToBinary(array)}
    winners = []
    binaries.each_with_index do |num, i|
      if num == binaries.min then winners.push(pcSetArray[i]) end
    end
    winners.sort[0]
  end
  
end