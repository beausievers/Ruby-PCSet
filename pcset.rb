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

def array_to_binary(array)
  array.inject(0) {|sum, n| sum + 2**n}
end

# the following method is horrifically inelegant 
# but avoids monkey-patching.
# TODO: do this right, incl. error checking
def pearsons(x, y)
  if !x.is_a?(Array) || !y.is_a?(Array) then raise StandardError, "x and y must be arrays", caller end
  if x.size != y.size then raise StandardError, "x and y must be same size", caller end
  sum_x = x.inject(0) {|sum, n| sum + n}
  sum_y = y.inject(0) {|sum, n| sum + n}
  sum_square_x = x.inject(0) {|sum, n| sum + n * n}
  sum_square_y = y.inject(0) {|sum, n| sum + n * n}
  xy = []
  x.zip(y) {|a, b| xy.push(a * b)}
  sum_xy = xy.inject(0) {|sum, n| sum + n}
  num = sum_xy - ((sum_x * sum_y)/x.size)
  den = Math.sqrt((sum_square_x - ((sum_x*sum_x)/x.size)) * (sum_square_y - ((sum_y*sum_y)/x.size)))
  (num/den)
end

class PCSet
  include Comparable
  attr_reader :pitches, :base, :input
  
  def initialize(pcarray, base = 12)
    if pcarray.instance_of?(Array) && pcarray.all?{|pc| pc.instance_of?(Fixnum)}
      @base, @input = base, pcarray
      @pitches = pcarray.map{ |x| x % @base }.uniq
    else
      raise ArgumentError, "Improperly formatted PC array", caller
    end
  end
  
  def PCSet.new_from_string(pcstring, base = 12)
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
  
  def invert(axis = 0)
    PCSet.new @pitches.map {|x| (axis-x) % @base}
  end
  
  def invert!(axis = 0)
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
  
  def transpositions_and_inversions(axis = 0)
    transpositions + invert(axis).transpositions
  end
  
  #
  # Normal form after Straus. Morris and AthenaCL do this differently.
  #
  def normal_form
    tempar = @pitches.sort
    arar = []                 # [[1,4,7,8,10],[4,7,8,10,1], etc.] get each cyclic variation
    tempar.each {arar.push PCSet.new(tempar.unshift(tempar.pop))}
    most_left_compact(arar)
  end
  
  def normal_form!
    @pitches = normal_form.pitches
  end
  
  def is_normal_form?
    self.pitches == self.normal_form.pitches
  end
  
  def set_class
    transpositions_and_inversions.map{|pcs| pcs.normal_form}.sort
  end
  
  def prime
    most_left_compact([normal_form.zero, invert.normal_form.zero])
  end
  
  def prime!
    self.pitches = self.prime.pitches
  end
  
  def is_prime?
    self.pitches == self.prime.pitches
  end
  
  def complement
    new_pitches = []
    @base.times do |p|
      if !@pitches.include? p then
        new_pitches.push p
      end
    end
    PCSet.new new_pitches
  end
  
  def full_interval_vector
    pairs = choose(@pitches, 2)                       # choose every pc pair
    intervals = pairs.map {|x| (x[1] - x[0]) % @base} # calculate every interval
    i_vector = Array.new(@base-1).fill(0)
    intervals.each {|x| i_vector[x-1] += 1}            # count the intervals
    i_vector
  end
  
  def interval_vector
    i_vector = full_interval_vector
    (0..((@base-1)/2)-1).each {|x| i_vector[x] += i_vector.pop}
    i_vector
  end
  
  #
  # Morris's invariance vector
  #
  def invariance_vector(m = 5)
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
    huron_table = [-1.428, -0.582, 0.594, 0.386, 1.240, -0.453]
    interval_consonance = []
    interval_vector.zip(huron_table) {|x, y| interval_consonance.push(x * y) }
    aggregate_dyadic_consonance = interval_consonance.inject {|sum, n| sum + n}
    [aggregate_dyadic_consonance, pearsons(interval_vector, huron_table)]
  end

  #
  # Balzano's vector of relations. Citation for all Balzano methods:
  # 
  # Balzano. "The Pitch Set as a Level of Description for Studying Musical
  # Pitch Perception" in Music, Mind, and Brain ed. Clynes. Plenum Press. 1982.
  #
  def vector_of_relations
    (0..length-1).to_a.map do |i|
      (0..length-1).to_a.map do |j|
        (@pitches[(i + j) % length] - @pitches[i]) % @base
      end
    end
  end

  #
  # Checks if the set satisfies Balzano's uniqueness.
  #
  def is_unique?
    vector_of_relations.uniq.size == vector_of_relations.size
  end
  
  #
  # Checks if the set satisfies Balzano's scalestep-semitone coherence.
  # For all s[i] and s[i1]:
  # j < k => v[i][j] < v[i1][k]
  # Where j and k are scalestep-counting indices.
  # And unless v[i][j] == 6 (a tritone), in which case the strict inequality is relaxed.
  #
  def is_coherent?
    v = vector_of_relations
    truth_array = []
    all_pair_indices = choose((0..length-1).to_a, 2)
    all_pair_indices.each do |i, i1|
      all_pair_indices.each do |j, k|
        if v[i][j] == 6
          truth_array.push(v[i][j] <= v[i1][k])
        else
          truth_array.push(v[i][j] < v[i1][k])    
        end
        if v[i1][j] == 6
          truth_array.push(v[i1][j] <= v[i][k])
        else
          truth_array.push(v[i1][j] < v[i][k])
        end
      end
    end
    !truth_array.include?(false)
  end
  
  #
  # Strict Balzano coherence, no inequality relaxation for tritones.
  #
  def is_strictly_coherent?
    v = vector_of_relations
    truth_array = []
    all_pair_indices = choose((0..length-1).to_a, 2)
    all_pair_indices.each do |i, i1|
      all_pair_indices.each do |j, k|
        truth_array.push(v[i][j] < v[i1][k])
        truth_array.push(v[i1][j] < v[i][k])
      end
    end
    !truth_array.include?(false)
  end
  
  def notes(middle_c = 0)
    noteArray = ['C','C#','D','D#','E','F','F#','G','G#','A','A#','B']
    if @base != 12 then raise StandardError, "PCSet.notes only makes sense for mod 12 pcsets", caller end
    out_string = String.new
    transpose(-middle_c).pitches.each do |p|
      out_string += noteArray[p] + ", "
    end
    out_string.chop.chop
  end

  def info
    print "modulo: #{@base}\n"
    print "raw input: #{@input.inspect}\n"
    print "pitch set: #{@pitches.inspect}\n"
    print "notes: #{notes}\n"
    print "normal: #{normal_form.inspect}\n"
    print "prime: #{prime.inspect}\n"
    print "interval vector: #{interval_vector.inspect}\n"
    print "invariance vector: #{invariance_vector.inspect}\n"
    print "huron ADC: #{huron[0]}  pearsons: #{huron[1]}\n"
    print "balzano coherence: "
    if is_strictly_coherent?
      print "strictly coherent\n"
    elsif is_coherent?
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
  def most_left_compact(pcset_array)
    if !pcset_array.all? {|pcs| pcs.length == pcset_array[0].length}
      raise ArgumentError, "PCSet.most_left_compact: All PCSets must be of same cardinality", caller
    end
    zeroed_pitch_arrays = pcset_array.map {|pcs| pcs.zero.pitches}
    binaries = zeroed_pitch_arrays.map {|array| array_to_binary(array)}
    winners = []
    binaries.each_with_index do |num, i|
      if num == binaries.min then winners.push(pcset_array[i]) end
    end
    winners.sort[0]
  end
  
end