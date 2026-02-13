@0xcce15f4779921ec4; #file id 

using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("panman");

struct NucMut
{
    nucPosition @0: Int32;
    nucGapPosition @1: Int32;  
    nucGapExist @2: Bool;
    mutInfo @3: UInt32; 
}

struct Mutation
{
    chrIdx @6: Int64;
    blockId @0: Int64;
    blockGapExist @1: Bool;
    blockMutExist @2: Bool;
    blockMutInfo @3: Bool;
    blockInversion @4: Bool;
    nucMutation @5: List(NucMut);
}

struct Node
{
    mutations @0: List(Mutation);
    annotations @1: List(Text);
}

struct ConsensusSeqToBlockIds
{
    blockId @0: List(Int64);
    consensusSeq @1: List(UInt32);
    blockGapExist @2: List(Bool);
    chromosomeName @3: List(Text);
}

struct GapList
{
    blockId @0: Int64;
    blockGapExist @1: Bool;
    nucGapLength @2: List(Int32);
    nucPosition @3: List(Int32);
}

struct BlockGapList
{
    blockPosition @0: List(Int32);
    blockGapLength @1: List(Int32);
}

struct ChrList
{
    chrIdx @0: Int64;
    chrName @1: Text;
    blockIds @2: List(Int64);
}

struct CircularOffset
{
    sequenceId @0: Text;
    offset @1: Int32;
}

struct RotationIndex
{
    sequenceId @0: Text;
    blockOffset @1: Int32;
}

struct SequenceInverted
{
    sequenceId @0: Text;
    inverted @1: Bool;
}

struct Tree
{
    newick @0: Text;
    nodes @1: List(Node);
    consensusSeqMap @2: List(ConsensusSeqToBlockIds);
    gaps @3: List(GapList);
    blockGaps @4: BlockGapList;
    circularSequences @5: List(CircularOffset);
    rotationIndexes @6: List(RotationIndex);
    sequencesInverted @7: List(SequenceInverted);
    chrLists @8: List(ChrList);
}

struct ComplexMutation {
    mutationType @0: Bool;
    treeIndex1 @1: Int32;
    treeIndex2 @2: Int32;
    treeIndex3 @3: Int32;
    sequenceId1 @4: Text;
    sequenceId2 @5: Text;

    blockIdStart1 @6: Int64;
    blockGapExistStart1 @7: Bool;
    nucPositionStart1 @8: Int32;
    nucGapPositionStart1 @9: Int32;
    nucGapExistStart1 @10: Bool;

    blockIdEnd1 @11: Int64;
    blockGapExistEnd1 @12: Bool;
    nucPositionEnd1 @13: Int32;
    nucGapPositionEnd1 @14: Int32;
    nucGapExistEnd1 @15: Bool;

    blockIdStart2 @16: Int64;
    blockGapExistStart2 @17: Bool;
    nucPositionStart2 @18: Int32;
    nucGapPositionStart2 @19: Int32;
    nucGapExistStart2 @20: Bool;

    blockIdEnd2 @21: Int64;
    blockGapExistEnd2 @22: Bool;
    nucPositionEnd2 @23: Int32;
    nucGapPositionEnd2 @24: Int32;
    nucGapExistEnd2 @25: Bool;

    sequenceId3 @26: Text;
}

struct TreeGroup
{
    trees @0: List(Tree);
    complexMutations @1: List(ComplexMutation);
}
