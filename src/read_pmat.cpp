#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <google/protobuf/message.h>

#include "PangenomeMAT.hpp"

#include "mutation_annotation_test_proto3_optional.pb.h"
#include "panman.pb.h"
bool compare_pair(const std::pair<std::pair<int,int>, int> &a, const std::pair<std::pair<int,int>, int> &b)
{
    return (a.first.first < b.first.first);
}

bool compare_pair_pair(const std::pair<std::pair<int,int>, int> &a, std::pair<std::pair<int,int>, int> &b)
{
    if (a.first.first != b.first.first)
        return (a.first < b.first);
    
    return (a.first.second < b.first.second);
}

/*
void old_proto_to_new_proto(MATNew::tree &tree)
{
    panman::tree panmanTree;

    panmanTree.set_newick(tree.newick());
    for (auto &n: tree.nodes())
    {
        std::vector<std::pair<std::pair<int,int>, int>> insert_blockid_global;
        std::vector<std::pair<std::pair<int,int>, int>> insert_blockid_local;
        std::vector<std::pair<std::pair<int,int>, int>> del_blockid_global;
        std::vector<std::pair<std::pair<int,int>, int>> del_blockid_local;
        std::vector<std::pair<std::pair<int,int>, int>> noindel_blockid_global;
        std::vector<std::pair<std::pair<int,int>, int>> noindel_blockid_local;
        std::unordered_map<int, std::vector<MATNew::nucMut>> nuc_mutations_list;
        std::unordered_map<int, bool> inversion_list;

        
        
        int idx = 0;
        for (auto &bid: n.mutations())
        {

            int64_t blockid = bid.blockid();
            int32_t global_idx = blockid>>32;
            int32_t local_idx = blockid;
            inversion_list[idx] = bid.blockinversion();
            if (!bid.blockmutexist())
            {
                if (bid.blockgapexist())
                    noindel_blockid_global.push_back(std::make_pair(std::make_pair(global_idx, local_idx), idx));
                else
                    noindel_blockid_local.push_back(std::make_pair(std::make_pair(global_idx, local_idx), idx));
            }
            else if (bid.blockmutinfo())
            {
                if (bid.blockgapexist())
                    insert_blockid_global.push_back(std::make_pair(std::make_pair(global_idx, local_idx), idx));
                else
                    insert_blockid_local.push_back(std::make_pair(std::make_pair(global_idx, local_idx), idx));
            }
            else
            {
                if (bid.blockgapexist())
                    del_blockid_global.push_back(std::make_pair(std::make_pair(global_idx, local_idx), idx));
                else
                    del_blockid_local.push_back(std::make_pair(std::make_pair(global_idx, local_idx), idx));
            }

            for (auto &nucm: bid.nucmutation())
            {
                nuc_mutations_list[idx].push_back(nucm); 
                
            }
            idx ++;
        }

        // sorting global co-ordinates
        std::sort(insert_blockid_global.begin(), insert_blockid_global.end(),compare_pair);
        std::sort(del_blockid_global.begin(), del_blockid_global.end(),compare_pair);
        std::sort(noindel_blockid_global.begin(), noindel_blockid_global.end(),compare_pair);

        // sorting local co-ordinates
        std::sort(insert_blockid_local.begin(), insert_blockid_local.end(), compare_pair_pair);
        std::sort(del_blockid_local.begin(), del_blockid_local.end(),compare_pair_pair);
        std::sort(noindel_blockid_local.begin(), noindel_blockid_local.end(),compare_pair_pair);


        panman::node newNode;

        for (auto &bid: n.annotations())
        {
            newNode.add_annotations();
            *newNode.mutable_annotations(newNode.annotations_size() - 1) = bid;
        }


        int length = 0;
        int id = 0;
        int current_id = 0;
        int64_t block_id;

        for (auto i = 0; i < noindel_blockid_global.size() + 1; i++)
        {
            if (i < noindel_blockid_global.size())
            {
                if ( length == 0 )
                {
                    block_id = noindel_blockid_global[0].first.first<<32 | noindel_blockid_global[0].first.second;
                    id = noindel_blockid_global[0].first.first;
                    length++;
                    continue;
                }
                current_id = noindel_blockid_global[0].first.first;
                if (current_id == id + 1)
                {
                    length++;
                    id = current_id;
                    continue;
                }
            }

            int8_t value = 0;
            value |= (length << 3);
            
            char charValue[0];
            charValue[0] = static_cast<char>(value);
            const char* char_pointer = charValue; 

            panman::mutation mutate;
            mutate.set_blockid(block_id);
            mutate.set_blockmutinfo(char_pointer);

            for (auto j = length - 1; j >= 0; j--)
            {
                panman::nucMut nucMut;
                int map_index = noindel_blockid_global[i - j - 1].second;

                for (auto &z: nuc_mutations_list[map_index])
                {
                    panman::nucMutLocal nucMutLocal;
                    nucMutLocal.set_nucposition(z.nucposition());
                    nucMutLocal.set_nucgapposition(z.nucgapposition());
                    nucMutLocal.set_nucgapexist(z.nucgapexist());
                    nucMutLocal.set_mutinfo(z.mutinfo());

                    nucMut.add_nucmutation();
                    *nucMut.mutable_nucmutation(nucMut.nucmutation_size() - 1) = nucMutLocal;
                }
                nucMut.set_blockinversion(inversion_list[map_index]);

                mutate.add_nucmutations();
                *mutate.mutable_nucmutations(mutate.nucmutations_size() - 1) = nucMut;
            }

            length = 1;
            if (i < noindel_blockid_global.size())
            {
                block_id = noindel_blockid_global[i].first.first<<32 | noindel_blockid_global[i].first.second;
                id = noindel_blockid_global[i].first.first;
            }

            
            newNode.add_mutations();
            *newNode.mutable_mutations(newNode.mutations_size() - 1) = mutate;
        } 

        
        length = 0;
        id = 0;
        current_id = 0;
        for (auto i = 0; i < insert_blockid_global.size() + 1; i++)
        {
            if (i < insert_blockid_global.size())
            {
                if ( length == 0 )
                {
                    block_id = insert_blockid_global[0].first.first<<32 | insert_blockid_global[0].first.second;
                    id = insert_blockid_global[0].first.first;
                    length++;
                    continue;
                }
                current_id = insert_blockid_global[0].first.first;
                if (current_id == id + 1)
                {
                    length++;
                    id = current_id;
                    continue;
                }
            }
            int8_t value = 0;
            value |= (length << 3);
            value |= (1 << 1); //mut-exist
            value |= (1 << 0); //mut-info -> 1 -> Insertion

            char charValue[0];
            charValue[0] = static_cast<char>(value);
            const char* char_pointer = charValue; 

            panman::mutation mutate;
            mutate.set_blockid(block_id);
            mutate.set_blockmutinfo(char_pointer);

            for (auto j = length - 1; j >= 0; j--)
            {
                panman::nucMut nucMut;
                int map_index = insert_blockid_global[i - j - 1].second;

                for (auto &z: nuc_mutations_list[map_index])
                {
                    panman::nucMutLocal nucMutLocal;
                    nucMutLocal.set_nucposition(z.nucposition());
                    nucMutLocal.set_nucgapposition(z.nucgapposition());
                    nucMutLocal.set_nucgapexist(z.nucgapexist());
                    nucMutLocal.set_mutinfo(z.mutinfo());

                    nucMut.add_nucmutation();
                    *nucMut.mutable_nucmutation(nucMut.nucmutation_size() - 1) = nucMutLocal;
                }
                nucMut.set_blockinversion(inversion_list[map_index]);

                mutate.add_nucmutations();
                *mutate.mutable_nucmutations(mutate.nucmutations_size() - 1) = nucMut;
            }

            length = 1;
            if (i < insert_blockid_global.size())
            {
                block_id = insert_blockid_global[i].first.first<<32 | insert_blockid_global[i].first.second;
                id = insert_blockid_global[i].first.first;

            }

            
            newNode.add_mutations();
            *newNode.mutable_mutations(newNode.mutations_size() - 1) = mutate;
        }     

        length = 0;
        id = 0;
        current_id = 0;
        for (auto i = 0; i < del_blockid_global.size() + 1; i++)
        {
            if (i < del_blockid_global.size())
            {
                if ( length == 0 )
                {
                    block_id = del_blockid_global[0].first.first<<32 | del_blockid_global[0].first.second;
                    id = del_blockid_global[0].first.first;
                    length++;
                    continue;
                }
                current_id = del_blockid_global[0].first.first;
                if (current_id == id + 1)
                {
                    length++;
                    id = current_id;
                    continue;
                }
            }

            int8_t value = 0;
            value |= (length << 3);
            value |= (1 << 1); //mut-exist
            value |= (0 << 0); //mut-info -> 0 -> Insertion

            char charValue[0];
            charValue[0] = static_cast<char>(value);
            const char* char_pointer = charValue; 

            panman::mutation mutate;
            mutate.set_blockid(block_id);
            mutate.set_blockmutinfo(char_pointer);

            for (auto j = length - 1; j >= 0; j--)
            {
                panman::nucMut nucMut;
                int map_index = del_blockid_global[i - j - 1].second;

                for (auto &z: nuc_mutations_list[map_index])
                {
                    panman::nucMutLocal nucMutLocal;
                    nucMutLocal.set_nucposition(z.nucposition());
                    nucMutLocal.set_nucgapposition(z.nucgapposition());
                    nucMutLocal.set_nucgapexist(z.nucgapexist());
                    nucMutLocal.set_mutinfo(z.mutinfo());

                    nucMut.add_nucmutation();
                    *nucMut.mutable_nucmutation(nucMut.nucmutation_size() - 1) = nucMutLocal;
                }
                nucMut.set_blockinversion(inversion_list[map_index]);
                mutate.add_nucmutations();
                *mutate.mutable_nucmutations(mutate.nucmutations_size() - 1) = nucMut;
            }

            length = 1;
            if (i < del_blockid_global.size())
            {
                block_id = del_blockid_global[i].first.first<<32 | del_blockid_global[i].first.second;
                id = del_blockid_global[i].first.first;
            }

            
            newNode.add_mutations();
            *newNode.mutable_mutations(newNode.mutations_size() - 1) = mutate;
        }   

        
        length = 0;
        std::pair<int,int> id_pair;
        std::pair<int,int> current_id_pair;
        for (auto i = 0; i < noindel_blockid_local.size() + 1; i++)
        {
            if (i < noindel_blockid_local.size())
            {
                if ( length == 0 )
                {
                    block_id = noindel_blockid_local[0].first.first<<32 | noindel_blockid_local[0].first.second;
                    id_pair = noindel_blockid_local[0].first;
                    length++;
                    continue;
                }
                current_id_pair = noindel_blockid_local[0].first;
                if ((current_id_pair.first == id_pair.first) && (current_id_pair.second == id_pair.second + 1))
                {
                    id_pair = current_id_pair;
                    length++;
                    continue;
                }
            }
            int8_t value = 0;
            value |= (length << 3);
            value |= (1 << 2); // gap_exit

            char charValue[0];
            charValue[0] = static_cast<char>(value);
            const char* char_pointer = charValue; 

            panman::mutation mutate;
            mutate.set_blockid(block_id);
            mutate.set_blockmutinfo(char_pointer);

            for (auto j = length - 1; j >= 0; j--)
            {
                panman::nucMut nucMut;
                int map_index = noindel_blockid_local[i - j - 1].second;

                for (auto &z: nuc_mutations_list[map_index])
                {
                    panman::nucMutLocal nucMutLocal;
                    nucMutLocal.set_nucposition(z.nucposition());
                    nucMutLocal.set_nucgapposition(z.nucgapposition());
                    nucMutLocal.set_nucgapexist(z.nucgapexist());
                    nucMutLocal.set_mutinfo(z.mutinfo());

                    nucMut.add_nucmutation();
                    *nucMut.mutable_nucmutation(nucMut.nucmutation_size() - 1) = nucMutLocal;
                }
                nucMut.set_blockinversion(inversion_list[map_index]);
                mutate.add_nucmutations();
                *mutate.mutable_nucmutations(mutate.nucmutations_size() - 1) = nucMut;
            }

            length = 1;
            if (i < noindel_blockid_local.size())
            {
                block_id = noindel_blockid_local[i].first.first<<32 | noindel_blockid_local[i].first.second;
                id_pair = noindel_blockid_local[i].first;
            }

            
            newNode.add_mutations();
            *newNode.mutable_mutations(newNode.mutations_size() - 1) = mutate;
        } 

        
        length = 0;
        for (auto i = 0; i < insert_blockid_local.size() + 1; i++)
        {
            if (i < insert_blockid_local.size())
            {
                if ( length == 0 )
                {
                    block_id = insert_blockid_local[0].first.first<<32 | insert_blockid_local[0].first.second;
                    id_pair = insert_blockid_local[0].first;
                    length++;
                    continue;
                }
                current_id_pair = insert_blockid_local[0].first;
                if ((current_id_pair.first == id_pair.first) && (current_id_pair.second == id_pair.second + 1))
                {
                    id_pair = current_id_pair;
                    length++;
                    continue;
                }
            }

            int8_t value = 0;
            value |= (length << 3);
            value |= (1 << 2); // gap_exit
            value |= (1 << 1); //mut-exist
            value |= (1 << 0); //mut-info -> 1 -> Insertion

            char charValue[0];
            charValue[0] = static_cast<char>(value);
            const char* char_pointer = charValue; 

            panman::mutation mutate;
            mutate.set_blockid(block_id);
            mutate.set_blockmutinfo(char_pointer);

            for (auto j = length - 1; j >= 0; j--)
            {
                panman::nucMut nucMut;
                int map_index = insert_blockid_local[i - j - 1].second;

                for (auto &z: nuc_mutations_list[map_index])
                {
                    panman::nucMutLocal nucMutLocal;
                    nucMutLocal.set_nucposition(z.nucposition());
                    nucMutLocal.set_nucgapposition(z.nucgapposition());
                    nucMutLocal.set_nucgapexist(z.nucgapexist());
                    nucMutLocal.set_mutinfo(z.mutinfo());

                    nucMut.add_nucmutation();
                    *nucMut.mutable_nucmutation(nucMut.nucmutation_size() - 1) = nucMutLocal;
                }
                nucMut.set_blockinversion(inversion_list[map_index]);
                mutate.add_nucmutations();
                *mutate.mutable_nucmutations(mutate.nucmutations_size() - 1) = nucMut;
            }

            length = 1;
            if (i < insert_blockid_local.size())
            {
                block_id = insert_blockid_local[i].first.first<<32 | insert_blockid_local[i].first.second;
                id_pair = insert_blockid_local[i].first;
            }

            
            newNode.add_mutations();
            *newNode.mutable_mutations(newNode.mutations_size() - 1) = mutate;
        }     

        length = 0;
        for (auto i = 0; i < del_blockid_local.size() + 1; i++)
        {
            if (i < del_blockid_local.size())
            {
                if ( length == 0 )
                {
                    block_id = del_blockid_local[0].first.first<<32 | del_blockid_local[0].first.second;
                    id_pair = del_blockid_local[0].first;
                    length++;
                    continue;
                }
                current_id_pair = del_blockid_local[0].first;
                if ((current_id_pair.first == id_pair.first) && (current_id_pair.second == id_pair.second + 1))
                {
                    id_pair = current_id_pair;
                    length++;
                    continue;
                }
            }

            int8_t value = 0;
            value |= (length << 3);
            value |= (1 << 2); // gap_exit
            value |= (1 << 1); //mut-exist
            value |= (0 << 0); //mut-info -> 0 -> Insertion

            char charValue[0];
            charValue[0] = static_cast<char>(value);
            const char* char_pointer = charValue; 

            panman::mutation mutate;
            mutate.set_blockid(block_id);
            mutate.set_blockmutinfo(char_pointer);

            for (auto j = length - 1; j >= 0; j--)
            {
                panman::nucMut nucMut;
                int map_index = del_blockid_local[i - j - 1].second;

                for (auto &z: nuc_mutations_list[map_index])
                {
                    panman::nucMutLocal nucMutLocal;
                    nucMutLocal.set_nucposition(z.nucposition());
                    nucMutLocal.set_nucgapposition(z.nucgapposition());
                    nucMutLocal.set_nucgapexist(z.nucgapexist());
                    nucMutLocal.set_mutinfo(z.mutinfo());

                    nucMut.add_nucmutation();
                    *nucMut.mutable_nucmutation(nucMut.nucmutation_size() - 1) = nucMutLocal;
                }
                nucMut.set_blockinversion(inversion_list[map_index]);
                mutate.add_nucmutations();
                *mutate.mutable_nucmutations(mutate.nucmutations_size() - 1) = nucMut;
            }

            length = 1;
            if (i < del_blockid_local.size())
            {
                block_id = del_blockid_local[i].first.first<<32 | del_blockid_local[i].first.second;
                id_pair = del_blockid_local[i].first;
            }

            
            newNode.add_mutations();
            *newNode.mutable_mutations(newNode.mutations_size() - 1) = mutate;
        }

        panmanTree.add_nodes();
        *panmanTree.mutable_nodes(panmanTree.nodes_size() - 1) = newNode;

        
        // panman::node newNode;
        // for (auto &mut: n.mutations())
        // {
        //     panman::mutation mutate;
            
        //     int64_t value = mut.blockid();
        //     int64_t cvalue = (value << 32)|(value >> 32);
        //     mutate.set_blockid(value);
        //     mutate.set_blockgapexist(mut.blockgapexist());
        //     mutate.set_blockmutexist(mut.blockmutexist());
        //     mutate.set_blockmutinfo(mut.blockmutinfo());
        //     mutate.set_blockinversion(mut.blockinversion());

        //     for (auto &nucm: mut.nucmutation())
        //     {
        //         panman::nucMut nucMutation;
        //         nucMutation.set_mutinfo(nucm.mutinfo());
        //         nucMutation.set_nucposition(nucm.nucposition());
        //         nucMutation.set_nucgapexist(nucm.nucgapexist());
        //         nucMutation.set_nucgapposition(nucm.nucgapposition());

        //         mutate.add_nucmutation();
        //         *mutate.mutable_nucmutation(mutate.nucmutation_size() - 1) = nucMutation;
        //     }
        //     newNode.add_mutations();
        //     *newNode.mutable_mutations(newNode.mutations_size() - 1) = mutate;
        // }
        // panmanTree.add_nodes();
        // *panmanTree.mutable_nodes(panmanTree.nodes_size() - 1) = newNode;
        
    }

    for (auto &n: tree.consensusseqmap())
    {
        panman::consensusSeqToBlockIds cmap;
        for (auto bid: n.blockid())
        {
            cmap.add_blockid(bid);
        }
        for (auto cseq: n.consensusseq())
        {
            cmap.add_consensusseq(cseq);
        }
        for (auto bge: n.blockgapexist())
        {
            cmap.add_blockgapexist(bge);
        }
        panmanTree.add_consensusseqmap();
        *panmanTree.mutable_consensusseqmap(panmanTree.consensusseqmap_size() - 1) = cmap;
    }

    for (auto &g: tree.gaps())
    {
        panman::gapList gl;
        gl.set_blockid(g.blockid());
        gl.set_blockgapexist(g.blockgapexist());
        
        for (auto &ng: g.nucgaplength())
        {
            gl.add_nucgaplength(ng);
        }

        for (auto &np: g.nucposition())
        {
            gl.add_nucposition(np);
        }
        panmanTree.add_gaps();
        *panmanTree.mutable_gaps(panmanTree.gaps_size() - 1) = gl;
    }

    {
        MATNew::blockGapList bgl = tree.blockgaps();
        panman::blockGapList bgl_panman;
        for (auto &bp: bgl.blockposition())
        {
            bgl_panman.add_blockposition(bp);
        }
        for (auto &bl: bgl.blockgaplength())
        {
            bgl_panman.add_blockgaplength(bl);
        }
        *panmanTree.mutable_blockgaps() = bgl_panman;
    }

    for (auto &g: tree.circularsequences())
    {
        panman::circularOffset co;
        co.set_sequenceid(g.sequenceid());
        co.set_offset(g.offset());
        panmanTree.add_circularsequences();
        *panmanTree.mutable_circularsequences(panmanTree.circularsequences_size() - 1) = co;
    }

    std::ofstream fout("newmat1");
    panmanTree.SerializeToOstream(&fout);


}
*/

void new_proto_to_old_proto(panman::tree &tree)
{
    MATNew::tree old_panman;

    // set newick tree
    old_panman.set_newick(tree.newick());

    // set node mutations
    for (auto &n: tree.nodes())
    {
        // set annotations
        MATNew::node newNode;
        for (auto &bid: n.annotations())
        {
            newNode.add_annotations();
            *newNode.mutable_annotations(newNode.annotations_size() - 1) = bid;
        }
        
        // set mutations
        for (auto &bid: n.mutations())
        {
            int64_t blockid = bid.blockid();
            int32_t prim = blockid >> 32;
            int32_t seco = blockid;
            std::string  blockmutinfo = bid.blockmutinfo();
            const uint8_t* p = reinterpret_cast<const uint8_t*>(blockmutinfo.c_str());
            uint8_t value = p[0];
            int length = value>>3;
            bool  blockGapExist = (value >> 2) & 0x1;
            bool  blockMutExist = (value>>1) & 0x1; 
            bool  blockMutInfo = value & 0x1;

            int index = 0;
            for (auto &nucmutlist: bid.nucmutations())
            {
                MATNew::mutation mut;
                if (blockGapExist)
                {
                    mut.set_blockid((prim << 32) | (seco + index));
                }
                else
                {
                    mut.set_blockid(((prim + index) << 32) | seco);
                }
                mut.set_blockgapexist(blockGapExist); 
                mut.set_blockmutexist(blockMutExist); 
                mut.set_blockmutinfo(blockMutInfo);
                for (auto &z: nucmutlist.nucmutation())
                {
                    MATNew::nucMut nucmut_old;
                    nucmut_old.set_nucposition(z.nucposition());
                    nucmut_old.set_nucgapposition(z.nucgapposition());
                    nucmut_old.set_nucgapexist(z.nucgapexist());
                    nucmut_old.set_mutinfo(z.mutinfo());
                    mut.add_nucmutation();
                    *mut.mutable_nucmutation(mut.nucmutation_size() - 1) = nucmut_old;
                }
                index++;
                newNode.add_mutations();
                *newNode.mutable_mutations(newNode.mutations_size() - 1) = mut;
            }

        }

        old_panman.add_nodes();
        *old_panman.mutable_nodes(old_panman.nodes_size() - 1) = newNode;

    
    }


    // set consensus sequence map
    for (auto &n: tree.consensusseqmap())
    {
        MATNew::consensusSeqToBlockIds cmap;
        for (auto bid: n.blockid())
        {
            cmap.add_blockid(bid);
        }
        for (auto cseq: n.consensusseq())
        {
            cmap.add_consensusseq(cseq);
        }
        for (auto bge: n.blockgapexist())
        {
            cmap.add_blockgapexist(bge);
        }
        old_panman.add_consensusseqmap();
        *old_panman.mutable_consensusseqmap(old_panman.consensusseqmap_size() - 1) = cmap;
    }

    // set gap list
    for (auto &g: tree.gaps())
    {
        MATNew::gapList gl;
        gl.set_blockid(g.blockid());
        gl.set_blockgapexist(g.blockgapexist());
        
        for (auto &ng: g.nucgaplength())
        {
            gl.add_nucgaplength(ng);
        }

        for (auto &np: g.nucposition())
        {
            gl.add_nucposition(np);
        }
        old_panman.add_gaps();
        *old_panman.mutable_gaps(old_panman.gaps_size() - 1) = gl;
    }

    // set block gaps
    {
        panman::blockGapList bgl = tree.blockgaps();
        MATNew::blockGapList bgl_panman;
        for (auto &bp: bgl.blockposition())
        {
            bgl_panman.add_blockposition(bp);
        }
        for (auto &bl: bgl.blockgaplength())
        {
            bgl_panman.add_blockgaplength(bl);
        }
        *old_panman.mutable_blockgaps() = bgl_panman;
    }

    // set circular sequence information
    for (auto &g: tree.circularsequences())
    {
        MATNew::circularOffset co;
        co.set_sequenceid(g.sequenceid());
        co.set_offset(g.offset());
        old_panman.add_circularsequences();
        *old_panman.mutable_circularsequences(old_panman.circularsequences_size() - 1) = co;
    }

    std::ofstream fout("old_panman");
    old_panman.SerializeToOstream(&fout);
}



std::vector<int> combine_sorted_array(std::vector<int> &ins)
{
    if (ins.size()<=1)
        return ins;

    std::sort(ins.begin(), ins.end());
    std::vector<int> ans;

    for (auto &m: ins)
    {
        std::cout << m << " ";
    }
    std::cout << std::endl;

    int start = 0;
    int curr = ins[0];
    for (auto i = 1; i < ins.size(); i++)
    {
        if (ins[i] == curr + 1)
        {
            curr = ins[i];
            continue;
        }
        else
        {
            ans.push_back(i - start);
            start = i;
            curr = ins[i];
        }
    }

    return ans;
}
/*
void checkBlockMuts(MATNew::tree &tree)
{
    for (auto &n: tree.nodes())
    {
        std::vector<int> ins;
        std::vector<int> del;
        for (auto &m: n.mutations())
        { 
            if (m.blockmutexist())
            {
                if (m.blockmutinfo())
                    ins.push_back(m.blockid());
                else
                    del.push_back(m.blockid());
            }
        }
        
        std::vector<int> ins_count = combine_sorted_array(ins);
        std::vector<int> del_count = combine_sorted_array(del);
        std::cout << ins.size() << " " << del.size() << " ";
        std::cout << ins_count.size() << " " << del_count.size() << std::endl;
    }

}
*/

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " input_proto_file.pb" << std::endl;
        return 1;
    }

    // Read the protobuf file
    std::fstream input(argv[1], std::ios::in | std::ios::binary);
    if (!input) {
        std::cerr << "Failed to open input file." << std::endl;
        return 1;
    }

    // Create a message object of your Protobuf type
    MATNew::tree tree;
    tree.ParseFromIstream(&input); 
    // old_proto_to_new_proto(tree);
    
    // panman::tree tree;
    // tree.ParseFromIstream(&input);   
    // new_proto_to_old_proto(tree);

    // checkBlockMuts(tree);
    // exit(1);

    
    // Node Size
    int node_size = 0;
    int nucmut_size = 0;
    int blkmut_size = 0;
    for (auto &n: tree.nodes())
    {
        node_size += n.ByteSizeLong();
        for (auto &m: n.mutations())
        { 
            blkmut_size += m.ByteSizeLong();
            for (auto &c: m.nucmutation())
            {
                nucmut_size+=c.ByteSizeLong();
            } 
        }
    }
    std::cout << "Nodes Size (Accurate): " << node_size/(1024*1024) << " MB" << std::endl;
    std::cout << "NucMut Size (Accurate): " << nucmut_size/(1024*1024) << " MB" << std::endl;
    std::cout << "BlockMut Size (Accurate): " << (blkmut_size - nucmut_size)/(1024*1024) << " MB" << std::endl;


    // Consensus Size
    int consensus_size = 0;
    int block_id = 0;
    int cons_seq = 0;
    int nucgaplen = 0;
    int cname = 0;
    for (auto &n: tree.consensusseqmap())
    {
        consensus_size += n.ByteSizeLong();
        block_id += n.blockid().size()*8;
        cons_seq += n.consensusseq().size()*4;
        nucgaplen += n.blockgapexist().size();
        cname += n.chromosomename_size();
    }
    std::cout << "Consensus Map Size (Accurate): " << consensus_size/(1024*1024) << " MB" << std::endl;
    std::cout << "consensusseq Size (Accurate): " << cons_seq/(1024*1024) << " MB" << std::endl;
    std::cout << "BlockID size (Accurate): " << block_id/(1024) << " KB" << std::endl;
    std::cout << "BlockGapExist size (Accurate): " << nucgaplen/(1024) << " KB" << std::endl;
    std::cout << "CName size (Accurate): " << cname/(1024) << " KB" << std::endl;

    // Gaps Size
    int gap_size = 0;
    for (auto &n: tree.gaps())
    {
        gap_size += n.ByteSizeLong();
    }
    std::cout << "Gap Size (Accurate): " << gap_size/(1024*1024) << " MB" << std::endl;

    // Block Gap Size
    // std::cout << "Block Gap Size: " << tree.blockgaps().ByteSizeLong() << std::endl;
    

    return 0;
}
