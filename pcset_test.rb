require 'pcset.rb'
require 'test/unit'

#
# When possible, test cases are adapted from 
# Introduction to Post-Tonal Theory by Joseph N. Straus,
# unless obvious or otherwise noted.
#

class PCSetTest < Test::Unit::TestCase
  def test_init
    #assert_raise(ArgumentError) {PCSet.new []}
    assert_raise(ArgumentError) {PCSet.new [1, 2, 3, 'string']}
    assert_raise(ArgumentError) {PCSet.new "string"}
    assert_raise(ArgumentError) {PCSet.new [1, 2, 3.6, 4]}
    assert_equal([0, 1, 2, 9], PCSet.new([0, 1, 2, 33, 13]).pitches)
    assert_equal([3, 2, 1, 11, 10, 0], PCSet.newFromString('321bac').pitches)
    assert_equal([0,2,4,5,7,11,9], PCSet.new([12,2,4,5,7,11,9]).pitches)
    assert_nothing_raised() {PCSet.new []}
  end
  
  def test_inversion
    
  end
  
  def test_transposition
    
  end
  
  def test_multiplication
    
  end
  
  #
  # set                   normal                 prime                 forte #
  # 0,2,4,7,8,11          7,8,11,0,2,4           0,1,4,5,7,9           6-31
  # 0,1,2,4,5,7,11        11,0,1,2,4,5,7         0,1,2,3,5,6,8         7-Z36
  # 0,1,3,5,6,7,9,10,11   5,6,7,9,10,11,0,1,3    0,1,2,3,4,6,7,8,10    9-8
  #
  
  def test_normalForm
    testPC = PCSet.new [0,4,8,9,11]
    assert_kind_of(PCSet, testPC.normalForm)
    assert_equal([8,9,11,0,4], testPC.normalForm.pitches)
    assert_equal([10,1,4,6], PCSet.new([1,6,4,10]).normalForm.pitches)
    assert_equal([2,4,8,10], PCSet.new([10,8,4,2]).normalForm.pitches)
    assert_equal([7,8,11,0,2,4], PCSet.new([0,2,4,7,8,11]).normalForm.pitches)
    assert_equal([11,0,1,2,4,5,7], PCSet.new([0,1,2,4,5,7,11]).normalForm.pitches)
    assert_equal([5,6,7,9,10,11,0,1,3], PCSet.new([0,1,3,5,6,7,9,10,11]).normalForm.pitches)  
  end
  
  def test_primeForm
    assert_equal([0,1,2,6], PCSet.new([5,6,1,7]).prime.pitches)
    assert_equal([0,1,4], PCSet.new([2,5,6]).prime.pitches)
    assert_equal([0,1,4,5,7,9], PCSet.new([0,2,4,7,8,11]).prime.pitches)
    assert_equal([0,1,2,3,5,6,8], PCSet.new([0,1,2,4,5,7,11]).prime.pitches)
    assert_equal([0,1,2,3,4,6,7,8,10], PCSet.new([0,1,3,5,6,7,9,10,11]).prime.pitches)
  end
  
  def test_setClass
    testPcs = PCSet.new([2,5,6])
    testPrime = testPcs.prime
    assert_equal([
        [2,5,6], [3,6,7], [4,7,8], [5,8,9], [6,9,10], [7,10,11],
        [8,11,0],[9,0,1], [10,1,2],[11,2,3],[0,3,4],  [1,4,5],
        [6,7,10],[7,8,11],[8,9,0], [9,10,1],[10,11,2],[11,0,3],
        [0,1,4], [1,2,5], [2,3,6], [3,4,7], [4,5,8],  [5,6,9]
      ].sort, PCSet.new([2,5,6]).setClass.map{|x| x.pitches})
    assert_equal(testPcs.setClass.map{|x| x.pitches}, testPrime.setClass.map{|x| x.pitches})
  end
  
  def test_intervalVector
    assert_equal([2,1,2,1,0,0], PCSet.new([0,1,3,4]).intervalVector)
    assert_equal([2,5,4,3,6,1], PCSet.new([0,1,3,5,6,8,10]).intervalVector)
    assert_equal([0,6,0,6,0,3], PCSet.new([0,2,4,6,8,10]).intervalVector)
  end
  
  def test_complement
    assert_equal([6,7,8,9,10,11], PCSet.new([0,1,2,3,4,5]).complement.pitches)
    assert_equal([3,4,5], PCSet.new([0,1,2], 6).complement.pitches)
  end
  
  #
  # Test values from (Morris 1991), pages 105-111
  # Citation:
  # Morris. Class Notes for Atonal Music Theory
  # Lebanon, NH. Frog Peak Music, 1991.
  #
  def test_invarianceVector
    assert_equal([1,0,0,0,5,6,5,5],PCSet.new([0,2,5]).invarianceVector)
    assert_equal([2,2,2,2,6,6,6,6],PCSet.new([0,1,6,7]).invarianceVector)
    assert_equal([6,6,6,6,6,6,6,6],PCSet.new([0,2,4,6,8,10]).invarianceVector)
    assert_equal([1,0,0,0,0,0,0,0],PCSet.new([0,1,2,3,4,5,8]).invarianceVector)
    assert_equal([1,0,0,1,0,0,0,0],PCSet.new([0,1,2,3,5,6,8]).invarianceVector)
    assert_equal([12,12,12,12,0,0,0,0],PCSet.new([0,1,2,3,4,5,6,7,8,9,10,11]).invarianceVector)
  end
  
  #
  # Test values from (Huron 1994). Huron rounds, thus the 0.01 margin of error.
  # Citation:
  # Huron. Interval-Class Content in Equally Tempered Pitch-Class Sets: 
  # Common Scales Exhibit Optimum Tonal Consonance. 
  # Music Perception (1994) vol. 11 (3) pp. 289-305
  #
  def test_huron
    h1 = PCSet.new([0,1,2,3,4,5,6,7,8,9,10,11]).huron
    assert_in_delta(-0.2, h1[0], 0.01)
    assert_in_delta(0.21, h1[1], 0.01)
    h2 = PCSet.new([0,2,4,5,7,9,11]).huron
    assert_in_delta(4.76, h2[0], 0.01)
    assert_in_delta(0.62, h2[1], 0.01)
  end
  
  def test_coherence
    
  end
end