from treeswift import *
from sys import argv
from readJson import *
import time
import sys
import networkx as nx
import networkx.algorithms.dag as dag
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fjson = open(argv[1])
data = json.load(fjson)
fjson.close()

def normal_rotation(data):
    seqPaths, seqStrands = getPaths(data)
    blocks_len = getBlocksLen(data)

    seqList = []
    for seq in seqPaths:
        seqList.append(seq)

    s1 = seqPaths[seqList[0]]
    plot=False
    if plot:
        for z in range(1,len(s1)):
            s2 = seqPaths[seqList[z]]

            rows, cols = (len(s2), len(s1))
            matrix = np.zeros(shape=(rows, cols))

            x=np.array([])
            y=np.array([])

            for i in range (len(s1)):
                for j in range(len(s2)):
                    if (s1[i]==s2[j]):
                        # matrix[j][i] = 1
                        x = np.append(x,[i])
                        y = np.append(y,[len(s2)-j])


            plt.scatter(x,y)
            plt.savefig(str(z))
            plt.figure().clear()



    hist = True
    if hist:
        for z in range(10):
            s4 = seqPaths[seqList[z]]
            map_ ={}
            s=s4
            avg_blk_len = 0
            max_blk_len = 0
            for i in range(len(s)):
                if s[i] not in map_:
                    map_[s[i]] = 1
                else:
                    map_[s[i]] += 1
                if (blocks_len[s[i]] > max_blk_len):
                    max_blk_len = blocks_len[s[i]]
                avg_blk_len += blocks_len[s[i]]
            avg_blk_len /= len(blocks_len)

            p_x = []
            p_y = []
            for i in map_:
                if (map_[i]>2):
                    p_x.append(blocks_len[i])
                    p_y.append(map_[i])

            plt.title("Frequency of all repeated blocks")
            plt.xlabel("Block Size")
            plt.ylabel("Frequency")
        
            # plt.hist(p, 30)
            plt.scatter(p_x,p_y)

            text = "Block Sizes\nMin:"+str(min(p_x))+"\nAvg:"+str(int(avg_blk_len))+"\nMax:"+str(max_blk_len)
            plt.text(1800,10,text)
            # plt.plot([min(p_x),min(p_x)], [1,40])
            # min_text = "Min Block Length: "+str(min(p_x))
            # plt.text(min(p_x), 40, min_text)
            
            # plt.plot([avg_blk_len, avg_blk_len], [1,40])
            # avg_text = "Avg block len:"+str(int(avg_blk_len))
            # plt.text(avg_blk_len, 40, avg_text)

            # max_text = "Max Block Len:"+str(max_blk_len)
            # plt.text(2500, 20, max_text)

            save_file = "hist" + str(z)
            plt.savefig(save_file)
            plt.figure().clear()

def generate_k_mer(path, K):
    s1_k_mer = []
    
    for i in range(len(path) - K + 1):
        k_mer =[]
        for z in range(i,i+K):
            k_mer.append(path[z])
        s1_k_mer.append(k_mer)
    return s1_k_mer

def k_mer_based_rotation(data, K, M):
    seqPaths, seqStrands = getPaths(data)
    blocks_len = getBlocksLen(data)

    seqList = []
    for seq in seqPaths:
        seqList.append(seq)

    s1_k_mer = generate_k_mer(seqPaths[seqList[0]], K)
    s2_k_mer = generate_k_mer(seqPaths[seqList[2]], K)

    rows, cols = (len(s2_k_mer), len(s1_k_mer))
    matrix = np.zeros(shape=(rows, cols))

    x=np.array([])
    y=np.array([])

    for i in range(cols):
        for j in range(rows):
            count = 0
            for t in range(K):
                if (s1_k_mer[i][t]) in s2_k_mer[j]:
                    count += 1
            if (count >= M):
                x = np.append(x,[i])
                y = np.append(y,[rows-j])
    
    plt.scatter(x,y)
    plt.savefig("k-mer")
    plt.figure().clear()


k_mer_based_rotation(data,6,3)
    