#include "fitchSankoff.cuh"

char getNucleotideFromCode(int code) {
    switch(code) {
    case 1:
        return 'A';
    case 2:
        return 'C';
    case 4:
        return 'G';
    case 8:
        return 'T';
    case 5:
        return 'R';
    case 10:
        return 'Y';
    case 6:
        return 'S';
    case 9:
        return 'W';
    case 12:
        return 'K';
    case 3:
        return 'M';
    case 14:
        return 'B';
    case 13:
        return 'D';
    case 11:
        return 'H';
    case 7:
        return 'V';
    case 15:
        return 'N';
    default:
        return '-';
    }
}

void post_order_traversal(panmanUtils::Node* node, std::unordered_map<std::string, std::pair<int,int>>& node_id_map, int& id, int &leaf_id, int* leaf_or_not, int* child_map, int* order, int &max_order, std::unordered_map<int, std::string>& reverse_node_id_map, std::unordered_map<int, int>& internal_node_id_map, int& internal_node) {
    std::string node_name;
    int node_id;

    if (node->children.size() == 0) {
        node_name = node->identifier;
        node_id = id++;
        reverse_node_id_map[node_id] = node_name;
        node_id_map[node_name] = std::make_pair(node_id,leaf_id++);
        order[node_id] = 0;
        return;
    }
    
    for (auto child : node->children)
        post_order_traversal(child, node_id_map, id, leaf_id, leaf_or_not, child_map, order, max_order, reverse_node_id_map, internal_node_id_map, internal_node);
    

    node_name = node->identifier;
    node_id = id++;
    int internal_node_id = internal_node++;
    internal_node_id_map[node_id] = internal_node_id;
    reverse_node_id_map[node_id] = node_name;
    node_id_map[node_name] = std::make_pair(node_id, -1);
    bool first_child = false;
    int last_child_id = -1;
    int new_order = 0;
    for (auto child : node->children){
        std::string child_name = child->identifier;
        if (!first_child){
            first_child = true;
            leaf_or_not[node_id] = node_id_map[child_name].first;
        } else {
            child_map[last_child_id] = node_id_map[child_name].first;
        }
        last_child_id = node_id_map[child_name].first;
        if (new_order < order[last_child_id])
            new_order = order[last_child_id];
    }
    order[node_id] = new_order+1;
    if (new_order+1>max_order)
        max_order = new_order+1;

    return;
}

__device__
int char2int(char c)
{
    int s = 0;
    switch (c) {
        case 'A':
        return 1;
    case 'C':
        return 2;
    case 'G':
        return 4;
    case 'T':
        return 8;
    case 'R':
        return 5;
    case 'Y':
        return 10;
    case 'S':
        return 6;
    case 'W':
        return 9;
    case 'K':
        return 12;
    case 'M':
        return 3;
    case 'B':
        return 14;
    case 'D':
        return 13;
    case 'H':
        return 11;
    case 'V':
        return 7;
    case 'N':
        return 15;
    default:
        return 0;

    }
    return s;
}

__global__ 
void fs_fwd(
    char * d_seqs, 
    int * d_leaf_map, 
    int * d_leaf_or_not, 
    int * d_child_map, 
    int * d_order,
    int num_nodes,
    int num_leaves,
    int sites,
    int max_order,
    int32_t * d_ancestor[],
    int * d_internal_node_id_map)
{
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int bs = blockDim.x;
    int gs = gridDim.x;

    for (int i=bx;i<sites;i+=gs) //sites
    {
        int local_order = 0;
        // int ancestor_start_idx = 16*i*num_nodes;
        while(local_order <= max_order)
        {
            for (int j=tx; j<num_nodes; j+=bs) //nodes
            {
                int ancestor_start_idx2 = (d_internal_node_id_map[j])*16;
                int node_id = j;
                int order = d_order[node_id];
                if (order == local_order)
                {
                    int leaf_id = d_leaf_map[node_id];
                    int leaf_or_not = d_leaf_or_not[node_id];
                    // int v;
                    // if (leaf_or_not == -1) //leaf
                    // {
                    //     v = char2int(d_seqs[leaf_id*sites+i]);
                    //     for (int k=0;k<16;k++)
                    //     {
                    //         if (k==v)   d_ancestor[i][ancestor_start_idx2+k] = 0;
                    //         else        d_ancestor[i][ancestor_start_idx2+k] = INF;
                    //     }

                    // }
                    if (leaf_or_not>-1) //internal node
                    {
                        for (int k=0;k<16;k++)
                            d_ancestor[i][ancestor_start_idx2+k] = 0;
                        
                        bool found_min = false;
                        int child_id = leaf_or_not;
                        while (child_id != -1)
                        {
                            int child_start_idx = (d_internal_node_id_map[child_id])*16;
                            bool is_child_leaf = (d_leaf_or_not[child_id] == -1);
                            int leaf_v = -1;
                            if (is_child_leaf) //leaf
                                leaf_v=char2int(d_seqs[d_leaf_map[child_id]*sites+i]);
                            
                            for (int k=0;k<16;k++)
                            {
                                int min_value=INF;
                                if (is_child_leaf) 
                                    min_value = (k!=leaf_v);
                                else
                                {
                                    for (int l=0; l<16;l++) 
                                    {
                                        if (min_value > d_ancestor[i][child_start_idx+l]+(k!=l))
                                            min_value = d_ancestor[i][child_start_idx+l]+(k!=l);
                                    }
                                }
                                if (min_value < INF)
                                    d_ancestor[i][ancestor_start_idx2+k]+=min_value;
                                else 
                                    printf("Ideally should not happen.. Report to swalia@ucsd.edu, Thanks!! (Site: %d, Node: %d Child Node %d)\n", i, j, child_start_idx);
                            }
                            child_id = d_child_map[child_id];
                        }
                    }
                    // d_ancestor[ancestor_start_idx+node_id] = v;
                }
            }
            __syncthreads();
            local_order++;
        }
    }

}

__global__ 
void fs_bwd(
    int32_t *d_ancestor[],
    int8_t  *d_bwd_states,
    char * d_ref, 
    int * d_leaf_or_not, 
    int * d_child_map, 
    int * d_order,
    int num_nodes,
    int num_leaves,
    int sites,
    int max_order,
    char * d_seqs,
    int * d_leaf_map,
    int * d_internal_node_id_map)
{
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int bs = blockDim.x;
    int gs = gridDim.x;

    for (int i=bx;i<sites;i+=gs) //sites
    {
        int local_order = max_order;
        // int ancestor_start_idx = 16*i*num_nodes;
        int states_start_idx = i*num_nodes;
        while(local_order > -1)
        {
            for (int j=tx; j<num_nodes; j+=bs) //nodes
            {
                int node_id = j;
                int order = d_order[node_id];
                if (order == local_order)
                {
                    int leaf_or_not = d_leaf_or_not[node_id];
                    int v;
                    if (local_order == max_order) //root
                    {
                        v = char2int(d_ref[i]);
                        d_bwd_states[states_start_idx+j] = v;
                    }
                    if (leaf_or_not != -1) // internal nodes
                    {
                        v=d_bwd_states[states_start_idx+j];
                        int child_id = leaf_or_not;
                        while (child_id != -1)
                        {
                            int child_start_idx = (d_internal_node_id_map[child_id])*16;
                            int min_value = INF, min_ptr=-1;
                            
                            bool is_child_leaf = (d_leaf_or_not[child_id] == -1);
                            int leaf_v = -1;
                            if (is_child_leaf)
                                leaf_v = char2int(d_seqs[d_leaf_map[child_id]*sites+i]);
                            
                            for (int k=0;k<16;k++)
                            {
                                int value;
                                if (is_child_leaf) // leaf
                                    value = (k!=leaf_v);
                                else // internal node
                                    value = (k!=v)+d_ancestor[i][child_start_idx+k];
                                if (value<min_value)
                                {
                                    min_value = value;
                                    min_ptr = k;
                                }
                            }
                            d_bwd_states[states_start_idx+child_id] = min_ptr;
                            child_id = d_child_map[child_id];
                        }
                    }
                    // d_ancestor[ancestor_start_idx+node_id] = v;
                }
            }
            __syncthreads();
            local_order--;
        }
    }

}

__global__
void fs_assign_mut(
    int8_t * d_bwd_states,
    int * d_leaf_or_not, 
    int * d_child_map, 
    int * d_order,
    int num_nodes,
    int num_leaves,
    int sites,
    int max_order,
    int8_t * d_muts
){
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int bs = blockDim.x;
    int gs = gridDim.x;

    int type = -1;

    for (int i=bx;i<sites;i+=gs) //sites
    {
        int local_order = max_order;
        int states_start_idx = i*num_nodes;
        int mut_start_idx = i*num_nodes;
        while(local_order > -1)
        {
            for (int j=tx; j<num_nodes; j+=bs) //nodes
            {
                int node_id = j;
                int order = d_order[node_id];
                if (order == local_order)
                {
                    int leaf_or_not = d_leaf_or_not[node_id];
                    if(leaf_or_not != -1) // internal nodes
                    {
                        int32_t node_state = d_bwd_states[states_start_idx + j]; 
                        if (node_state==-1)
                            d_muts[mut_start_idx+j]=-1;
                        else {
                            int child_id = leaf_or_not;
                            while(child_id > -1)
                            {
                                int32_t child_state = d_bwd_states[states_start_idx+child_id];
                                if (node_state != child_state) 
                                {
                                    if (node_state == 0) // insertion
                                        type = 2;
                                    else if (child_state == 0) // deletion
                                        type = 1;
                                    else // subs
                                        type = 0;
                                    d_muts[mut_start_idx+child_id] = (type<<4 | child_state);
                                } 
                                else 
                                {
                                    d_muts[mut_start_idx+child_id] = -1;
                                }
                                child_id = d_child_map[child_id];
                            }
                        }
                        
                    }
                }
            }
            __syncthreads();
            local_order--;
        }
    }
}

void allocate_mem_and_run(std::unordered_map<std::string, std::string>& seqs, int* leaf_or_not, int* child_map, int* order, std::unordered_map<std::string, std::pair<int, int>>& node_id_map,utility::util* u, std::unordered_map<int, std::string>&reverse_node_id_map, std::unordered_map<int, int>& internal_node_id_map){
    FILE *file = fopen("mutations.txt", "w");
    std::string error;
    
    int num_seq = u->num_tips;
    int * h_leaf_map = new int[u->num_nodes];
    for (int i=0; i<u->num_nodes; i++) 
        h_leaf_map[i]=-1;

    for (auto s: node_id_map) 
    {
        int leaf_id             = s.second.second;
        int node_id             = s.second.first;
        if (leaf_id>-1)
            h_leaf_map[node_id] = leaf_id;
    }

    if (u->ref_name == ""){
        printf("Currently the program requires users to pass reference to perform fitch-Sankoff\n");
        exit(0);
    }
    
    int *d_internal_node_id_map, *h_internal_node_id_map;
    char *d_seqs,    *h_seqs,    *h_seqs_global;
    char *d_ref_seq, *h_ref_seq, *h_ref_seq_global;
    int32_t *d_ancestor[u->local_batch_size];
    int8_t  *d_bwd_states;
    // [u->local_batch_size];
    int8_t  *d_muts;
    // [u->local_batch_size];
    int *d_leaf_or_not = new int [u->num_nodes];
    int *d_child_map   = new int [u->num_nodes];
    int *d_order       = new int [u->num_nodes];
    int *d_leaf_map    = new int [u->num_nodes];
    
    // allocate memory on CPU
    // h_seqs_global    = (char *)malloc(u->global_batch_size*num_seq*sizeof(char));
    // h_ref_seq_global = (char *)malloc(u->global_batch_size*num_seq*sizeof(char));
    h_seqs           = (char *)malloc(u->local_batch_size*num_seq*sizeof(char));
    h_ref_seq        = (char *)malloc(u->local_batch_size*sizeof(char));
    h_internal_node_id_map = (int*) malloc(u->num_nodes*sizeof(int));

    for (auto n: internal_node_id_map)
        h_internal_node_id_map[n.first] = n.second;
    
    // allocate memory on GPU
    size_t freeMemory, totalMemory;
    cudaMemGetInfo(&freeMemory, &totalMemory);
    printf("============= Memory Usage on GPU ===============\n");
    printf("Total memory on GPU:%lf GB\n", (double)totalMemory/(1024*1024*1024));
    printf("Available memory: %lf GB\n", (double)freeMemory/(1024*1024*1024));
    
    printf("Allocating %lf GB for Seqs on GPU\n", (double)num_seq*u->local_batch_size/(1024*1024*1024));
    printf("Allocating %lf GB for Tree structure on GPU\n", (double)u->num_nodes*5*4/(1024*1024*1024));
    printf("Allocating %lf GB for Fitch forward states on GPU\n",((double)16*(u->num_nodes-u->num_tips)*u->local_batch_size*4 + u->local_batch_size * 8)/(1024*1024*1024));
    printf("Allocating %lf GB for Ref seq on GPU\n", (double)u->local_batch_size/(1024*1024*1024));
    printf("Allocating %lf GB for Fitch bwd states seq on GPU\n", (double)u->local_batch_size*u->num_nodes/(1024*1024*1024));
    printf("Allocating %lf GB for Fitch mutation seq on GPU\n", (double)u->num_nodes*u->local_batch_size/(1024*1024*1024)); 
    double total_usage = (double)num_seq*u->local_batch_size/(1024*1024*1024) + \
                         (double)u->num_nodes*5*4/(1024*1024*1024) + \
                         ((double)16*(u->num_nodes-u->num_tips)*u->local_batch_size*4 + u->local_batch_size *8)/(1024*1024*1024) + \
                         (double)u->local_batch_size/(1024*1024*1024) + \
                         (double)u->local_batch_size*u->num_nodes/(1024*1024*1024) + \
                         (double)u->num_nodes*u->local_batch_size/(1024*1024*1024);        

    printf("\nTotal Usage: %lf GB\n", total_usage);
    
    for (int i=0; i<u->local_batch_size; i++)
    {
        cudaMalloc(&d_ancestor[i], 16*(u->num_nodes-u->num_tips)*sizeof(int32_t));
        // cudaMalloc(&d_bwd_states[i], u->num_nodes*sizeof(int8_t));
        // cudaMalloc(&d_muts[i], u->num_nodes*sizeof(int8_t));
    }
    
    int32_t** d_vectorsArray;
    int8_t** d_bwd_states_array;
    int8_t** d_muts_array;
    cudaMalloc(&d_vectorsArray, u->local_batch_size * sizeof(int32_t*));
    // cudaMalloc(&d_bwd_states_array, u->local_batch_size * sizeof(int8_t*));
    // cudaMalloc(&d_muts_array, u->local_batch_size * sizeof(int8_t*));
    
    cudaMemcpy(d_vectorsArray, d_ancestor, u->local_batch_size * sizeof(int32_t*), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_bwd_states_array, d_bwd_states, u->local_batch_size * sizeof(int8_t*), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_muts_array, d_muts, u->local_batch_size * sizeof(int8_t*), cudaMemcpyHostToDevice);

    error = cudaGetErrorString(cudaGetLastError()); 
    if (error != "no error") 
    {
        printf("ERROR: cudaMallocPitch %s!\n", error.c_str());
        exit(0);
    }
    cudaMalloc(&d_seqs, num_seq*u->local_batch_size*sizeof(char));
    cudaMalloc(&d_leaf_or_not, u->num_nodes*sizeof(int));
    cudaMalloc(&d_child_map, u->num_nodes*sizeof(int));
    cudaMalloc(&d_order, u->num_nodes*sizeof(int));
    cudaMalloc(&d_leaf_map, u->num_nodes*sizeof(int));
    cudaMalloc(&d_ref_seq, u->local_batch_size*sizeof(char));
    cudaMalloc(&d_bwd_states, u->local_batch_size*u->num_nodes*sizeof(int8_t));
    cudaMalloc(&d_muts, u->num_nodes*u->local_batch_size*sizeof(int8_t));
    cudaMalloc(&d_internal_node_id_map, u->num_nodes*sizeof(int));

    cudaMemGetInfo(&freeMemory, &totalMemory);
    printf("Free memory after allocation: %lf GB\n", (double)freeMemory/(1024*1024*1024));
    printf("=================================================\n");

    // Requires one time transfer
    cudaMemcpy(d_leaf_or_not, leaf_or_not, u->num_nodes*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_child_map, child_map, u->num_nodes*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_order, order, u->num_nodes*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_leaf_map, h_leaf_map, u->num_nodes*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_internal_node_id_map, h_internal_node_id_map, u->num_nodes*sizeof(int), cudaMemcpyHostToDevice);
    
    error = cudaGetErrorString(cudaGetLastError()); 
    if (error != "no error") printf("ERROR: Cuda memcpy Trees structure %s!\n", error.c_str());
    

    size_t fitch_global_position = 0;
    while (fitch_global_position<u->ref_seq.size())
    {
        auto batch_start = std::chrono::high_resolution_clock::now();
        if (u->ref_seq.size()-fitch_global_position<u->global_batch_size)
            u->global_batch_size = u->ref_seq.size()-fitch_global_position;
        fprintf(stderr,"Handling data from position %ld to %d\n", fitch_global_position, fitch_global_position+u->global_batch_size);
        
        std::unordered_map<std::string, std::string> mutations;
        // Reading Sequences from Disk
        auto read_start = std::chrono::high_resolution_clock::now();
        read_seqs(u->seq_file_name, seqs, fitch_global_position, u->global_batch_size, 1);
        auto read_end = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds read_time = read_end - read_start;
        fprintf(stderr, "Reading Sequence took %lf mins\n", ((double)read_time.count()/1000000000)/60);
    
        size_t fitch_local_position = 0;
        while (fitch_local_position<u->global_batch_size)
        {
            if (fitch_local_position+u->local_batch_size>u->global_batch_size)
                u->local_batch_size = u->global_batch_size-fitch_local_position;
            fprintf(stderr, "%d...", fitch_local_position);
            for (auto s: node_id_map) {
                int leaf_id             = s.second.second;
                int node_id             = s.second.first;
                std::string node_name   = s.first;
                if (leaf_id>-1){
                    std::string seq = seqs[node_name];
                    for (int i=0;i<u->local_batch_size;i++)
                        h_seqs[leaf_id*u->local_batch_size+i]=seq[fitch_local_position+i];
                }
            }

            for (int i=0; i<u->local_batch_size;i++)
                h_ref_seq[i]=u->ref_seq[fitch_global_position+fitch_local_position+i];

            cudaMemcpy(d_seqs, h_seqs, num_seq*u->local_batch_size*sizeof(char), cudaMemcpyHostToDevice);
            cudaMemcpy(d_ref_seq, h_ref_seq, u->local_batch_size*sizeof(char), cudaMemcpyHostToDevice);
            
            error = cudaGetErrorString(cudaGetLastError()); 
            if (error != "no error") printf("ERROR: Cuda memcpy %s!\n", error.c_str());

            auto fwd_start = std::chrono::high_resolution_clock::now();
            fs_fwd<<<1024,1024>>>(d_seqs, d_leaf_map, d_leaf_or_not, d_child_map, d_order, u->num_nodes, seqs.size(), u->local_batch_size, u->max_order, d_vectorsArray, d_internal_node_id_map);
            cudaDeviceSynchronize();
            auto fwd_end = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds fwd_time = fwd_end - fwd_start;
            error = cudaGetErrorString(cudaGetLastError());
            if (error != "no error")
                printf("ERROR: After fs_fwd - %s!\n", error.c_str());

if (0)
{
    int16_t * h_ancestor = new int16_t[u->local_batch_size*u->num_nodes*16];    
    cudaMemcpy(h_ancestor,    d_ancestor,    (u->local_batch_size*u->num_nodes*16)*sizeof(int16_t), cudaMemcpyDeviceToHost);
    error = cudaGetErrorString(cudaGetLastError());
    if (error != "no error")
        printf("ERROR: After Copy to Host %s!\n", error.c_str());
    for(auto n:node_id_map){
        int node_id = n.second.first;
        printf("%s\n",n.first.c_str());
        for (int j=0;j<16;j++){
            printf("%d\t", h_ancestor[node_id*16+j]);
        }
        printf("\n");
    }
}
        

            auto bwd_start = std::chrono::high_resolution_clock::now();
            fs_bwd<<<1024,1024>>>(d_vectorsArray, d_bwd_states, d_ref_seq, d_leaf_or_not, d_child_map, d_order, u->num_nodes, seqs.size(), u->local_batch_size, u->max_order, d_seqs, d_leaf_map, d_internal_node_id_map);
            cudaDeviceSynchronize();
            auto bwd_end = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds bwd_time = bwd_end-bwd_start;

            error = cudaGetErrorString(cudaGetLastError());
            if (error != "no error")
                printf("ERROR: After fs_bwd - %s!\n", error.c_str());

if (0)
{
    int16_t * h_states = new int16_t[u->local_batch_size*u->num_nodes];    
    cudaMemcpy(h_states,    d_bwd_states,    (u->local_batch_size*u->num_nodes)*sizeof(int16_t), cudaMemcpyDeviceToHost);
    error = cudaGetErrorString(cudaGetLastError());
    if (error != "no error")
        printf("ERROR: After Copy to Host %s!\n", error.c_str());
    int position = 0;
    for(auto n:node_id_map){
        int node_id = n.second.first;
        int16_t v = h_states[node_id + position*u->num_nodes];
        printf("%d %s\n",node_id, n.first.c_str());
        printf("%d\t",v);
        
        printf("\n");
    }
}
           

            auto mut_start = std::chrono::high_resolution_clock::now();
            fs_assign_mut<<<1024,1024>>>(d_bwd_states, d_leaf_or_not, d_child_map, d_order, u->num_nodes, seqs.size(), u->local_batch_size, u->max_order, d_muts);
            cudaDeviceSynchronize();
            auto mut_end = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds mut_time = mut_end-mut_start;
            error = cudaGetErrorString(cudaGetLastError());
            if (error != "no error")
                printf("ERROR: After fs_assign_mut - %s!\n", error.c_str());

if(1)
{

    int8_t * h_muts = new int8_t[u->num_nodes*u->local_batch_size];
    cudaMemcpy(h_muts,    d_muts,    (u->num_nodes*u->local_batch_size)*sizeof(int8_t), cudaMemcpyDeviceToHost);

    error = cudaGetErrorString(cudaGetLastError());
    if (error != "no error")
        printf("ERROR: After Copy to Host %s!\n", error.c_str());
    for(auto n:node_id_map){
        int node_id = n.second.first;
        if (reverse_node_id_map[node_id] == "node_1")
            continue;
        std::string s = "";
        
        // s+= reverse_node_id_map[node_id];
        // s += ":\t";
        int count_mutations=0;
        for (int i=0; i<u->local_batch_size; i++)
        {
            int8_t v = h_muts[i*u->num_nodes + node_id];
            if (v!=-1)
            {
                count_mutations++;
                int type = v>>4;
                int8_t c = v & 0x0F;
                
                if (type == 0) s+= 'S';
                if (type == 1) s+= 'D';
                if (type == 2) s+= 'I';
                s += std::to_string(fitch_global_position+fitch_local_position+i);
                s += getNucleotideFromCode(c);
                s += '\t';
            }
        }
        // s += "\n";
        if (count_mutations>0)
        {
            // fprintf(file, "%s",s.c_str());
            if (mutations.find(reverse_node_id_map[node_id]) == mutations.end())
                mutations[reverse_node_id_map[node_id]]="";
            mutations[reverse_node_id_map[node_id]].append(s);
        }
    }
}
            fitch_local_position += u->local_batch_size;
        }
        fprintf(stderr,"%d...\n", fitch_local_position);
        auto batch_end = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds batch_time = batch_end - batch_start;
        fprintf(stderr,"Batch completed in %lf mins\n", ((double)batch_time.count()/1000000000)/60);

        for (int i=0; i<u->global_batch_size; i++)
        {
            char c = u->ref_seq[fitch_global_position+i];
            if (c=='-')
            {
                bool found = false;
                for (auto s: seqs) {
                    if (s.second[i] != '-')
                    {
                        u->consensus[fitch_global_position+i] = s.second[i];
                        found = true;
                        break;
                    }
                }
                if (found == false) 
                {
                    printf("Error: Position %d has all -\n", i);
                    exit(0);
                }
            }
        }

if (1)
{
    for (auto &s: mutations)
    {
        fprintf(file, "%s:\t", s.first.c_str());
        fprintf(file, "%s\n", s.second.c_str());
    }
}

        fitch_global_position += u->global_batch_size;
    }
    
if (1){
    // node_1 mutations
    std::string s = "node_1:\t";
    for (int i=0; i<u->ref_seq.size(); i++)
    {
        char child_state = u->ref_seq[i];
        char node_state = u->consensus[i];
        if (node_state != child_state) 
        {
            if (node_state == '-') // insertion
                s += 'I';
            else if (child_state == '-') // deletion
                s += 'D';
            else // subs
                s += 'S';
            s+=std::to_string(i);
            s+=child_state;
            s+='\t';
        } 
    }
    s += "\n";
    fprintf(file, "%s",s.c_str());
}
    fclose(file);

    file = fopen("consensus.txt", "w"); 
    fprintf(file, "%s",u->consensus.c_str());
    fclose(file);

    cudaFree(d_seqs);
    cudaFree(d_ref_seq);
    cudaFree(d_bwd_states);
    cudaFree(d_ancestor);
    cudaFree(d_vectorsArray);
    cudaFree(d_leaf_or_not);
    cudaFree(d_child_map);
    cudaFree(d_order);
    cudaFree(d_leaf_map);
    cudaFree(d_muts);
    // cudaFree(d_muts_array);
    // cudaFree(d_bwd_states_array);

    return;
}

void fitch_sankoff_on_gpu(panmanUtils::Tree* T, std::unordered_map<std::string, std::string>& seqs, utility::util* u) {
    printf("Creating Tree for Device\n");
    int *leaf_or_not,*child_map, *order;
    leaf_or_not = (int*)malloc(T->allNodes.size()*sizeof(int));
    child_map = (int*)malloc(T->allNodes.size()*sizeof(int));
    order = (int*)malloc(T->allNodes.size()*sizeof(int));
    for (int i = 0; i < T->allNodes.size(); i++) {
        child_map[i] = -1;
        leaf_or_not[i] = -1;
        order[i] = 0;
    }

    // create a mapping between node's old ids and new ids
    std::unordered_map<std::string, std::pair<int,int>> node_id_map;
    std::unordered_map<int, int> internal_node_id_map; // ids for nodes except leaf nodes
    std::unordered_map<int, std::string> reverse_node_id_map;
    int id = 0, leaf_id=0, max_order=0; int internal_node=0;
    // node_id_map[T->root->identifier.c_str()] = std::make_pair(id,-1);
    // internal_node_id_map[id] = internal_node;
    post_order_traversal(T->root, node_id_map, id, leaf_id, leaf_or_not, child_map, order, max_order, reverse_node_id_map, internal_node_id_map, internal_node);
    u->max_order = max_order;
    u->num_nodes = id;
    u->num_tips = leaf_id;
    
if (0) {
    for (auto node: T->allNodes) {
        int new_id = node_id_map[node.first.c_str()].first;
        printf("%d: %s, is a leaf: %d and has children: %d, and order: %d\n", new_id, node.first.c_str(), leaf_or_not[new_id], child_map[new_id], order[new_id]);
        printf("Position %d - char %c\n", 0, seqs[node.first][0]);
    }
}  

    allocate_mem_and_run(seqs, leaf_or_not, child_map, order, node_id_map,u, reverse_node_id_map, internal_node_id_map);

}
