/*
 * Copyright 1994-2011, Regents of the University of Minnesota
 *
 * gpmetis.c
 *
 * Drivers for the partitioning routines
 *
 * Started 8/28/94
 * George
 *
 * $Id: gpmetis.c 13900 2013-03-24 15:27:07Z karypis $
 *
 */

#include "metisbin.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

// A utility function to swap to integers
void swap_position (int *from, int *to)
{
    int temp = *from;
    *from = *to;
    *to = temp;
}

/*************************************************************************/
/*! Let the game begin! */
/*************************************************************************/
int main(int argc, char *argv[]) {
    printf("Method called\n");
    int argi=1;
    for(argi=1; argi < argc; ++argi){
        printf("[%d] param: %s \n", argi, argv[argi]);
    }
    idx_t i, j, k, l, cl, u, v, k_part;
    char *curptr, *newptr;
    idx_t options[METIS_NOPTIONS];
    graph_t *graph;
    idx_t * part;
    idx_t objval;
    params_t *params;
    int status = 0;

    params = parse_cmdline(argc, argv);

    gk_startcputimer(params->iotimer);
    graph = ReadGraph(params);
    printf("Read the matrix for %s\n", params->filename);
//    graph = ReadMatrix(params);
    ReadTPwgts(params, graph->ncon);
    gk_stopcputimer(params->iotimer);

    /* Check if the graph is contiguous */
    if (params->contig && !IsConnected(graph, 0)) {
        printf("***The input graph is not contiguous.\n"
               "***The specified -contig option will be ignored.\n");
        params->contig = 0;
    }

    /* Get ubvec if supplied */
    if (params->ubvecstr) {
        params->ubvec = rmalloc(graph->ncon, "main");
        curptr = params->ubvecstr;
        for (i = 0; i < graph->ncon; i++) {
            params->ubvec[i] = strtoreal(curptr, &newptr);
            if (curptr == newptr)
                errexit("Error parsing entry #%"PRIDX" of ubvec [%s] (possibly missing).\n",
                        i, params->ubvecstr);
            curptr = newptr;
        }
    }

    /* Setup iptype */
    if (params->iptype == -1) {
        if (params->ptype == METIS_PTYPE_RB) {
            if (graph->ncon == 1)
                params->iptype = METIS_IPTYPE_GROW;
            else
                params->iptype = METIS_IPTYPE_RANDOM;
        }
    }

    GPPrintInfo(params, graph);

    part = imalloc(graph->nvtxs, "main: part");

    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = params->objtype;
    options[METIS_OPTION_CTYPE] = params->ctype;
    options[METIS_OPTION_IPTYPE] = params->iptype;
    options[METIS_OPTION_RTYPE] = params->rtype;
    options[METIS_OPTION_NO2HOP] = params->no2hop;
    options[METIS_OPTION_MINCONN] = params->minconn;
    options[METIS_OPTION_CONTIG] = params->contig;
    options[METIS_OPTION_SEED] = params->seed;
    options[METIS_OPTION_NITER] = params->niter;
    options[METIS_OPTION_NCUTS] = params->ncuts;
    options[METIS_OPTION_UFACTOR] = params->ufactor;
    options[METIS_OPTION_DBGLVL] = params->dbglvl;

    gk_malloc_init();
    gk_startcputimer(params->parttimer);

    switch (params->ptype) {
        case METIS_PTYPE_RB:
            status = METIS_PartGraphRecursive(&graph->nvtxs, &graph->ncon, graph->xadj,
                                              graph->adjncy, graph->vwgt, graph->vsize, graph->adjwgt,
                                              &params->nparts, params->tpwgts, params->ubvec, options,
                                              &objval, part);
            break;

        case METIS_PTYPE_KWAY:
            status = METIS_PartGraphKway(&graph->nvtxs, &graph->ncon, graph->xadj,
                                         graph->adjncy, graph->vwgt, graph->vsize, graph->adjwgt,
                                         &params->nparts, params->tpwgts, params->ubvec, options,
                                         &objval, part);
            break;

    }

    gk_stopcputimer(params->parttimer);

    if (gk_GetCurMemoryUsed() != 0)
        printf("***It seems that Metis did not free all of its memory! Report this.\n");
    params->maxmemory = gk_GetMaxMemoryUsed();
    gk_malloc_cleanup(0);


    if (status != METIS_OK) {
        printf("\n***Metis returned with an error.\n");
    } else {
        if (!params->nooutput) {
            /* Write the solution */
            gk_startcputimer(params->iotimer);
            WritePartition(params->filename, part, graph->nvtxs, params->nparts);
            gk_stopcputimer(params->iotimer);
        }

        GPReportResults(params, graph, part, objval);
    }

    printf("Partition done!!!\n");
    //// Shuffle Graph ///
    /*idx_t *random_vartex = imalloc(graph->nvtxs, "main: part");
    for(u=0; u<graph->nvtxs; ++u){
        random_vartex[u] = u;
    }
    srand ( time(NULL) );
    for (i = graph->nvtxs-1; i > 0; i--) {
        int j = rand() % (i+1);
        // Swap arr[i] with the element at random index
        swap_position(&random_vartex[i], &random_vartex[j]);
    }

    idx_t new_id = 0, itr = 0;
    idx_t * new_ids;
    new_ids = imalloc(graph->nvtxs, "main: part");
    for (i = 0; i < graph->nvtxs; ++i) {
        new_ids[random_vartex[i]] = new_id++;
    }
    FILE *shuffleMat;
    char *ptr = strtok(params->filename, ".");

    char mat_filename[MAXLINE];
    sprintf(mat_filename, "%s_rm_shuffle", ptr);
    if (!(shuffleMat = fopen(strcat(mat_filename, ".mtx"), "w"))) {
        fprintf(stderr, "fopen: failed to open file '%s'", ptr);
        exit(EXIT_FAILURE);
    }
    fprintf(shuffleMat, "%%%MatrixMarket matrix coordinate real general\n");
    fprintf(shuffleMat, "%d %d %d\n", graph->nvtxs, graph->nvtxs, graph->nedges);

    for (i = 0; i < graph->nvtxs; ++i) {
        u = random_vartex[i];
        for (v = graph->xadj[u]; v < graph->xadj[u + 1]; v++) {
            fprintf(shuffleMat, "%d %d %lf\n", (new_ids[u] + 1), (new_ids[graph->adjncy[v]] + 1), (double) graph->adjwgt[v]);
        }
    }

    if (fclose(shuffleMat) != 0) {
        printf("fopen: failed to open file '%s'", ptr);
        exit(EXIT_FAILURE);
    }*/
    /// end shuffle /////

    /***** Randomize Matrix ******/
    idx_t *random_vartex = imalloc(graph->nvtxs, "main: part");
    idx_t row = sqrt(params->nparts);
    idx_t col = row;
    for(u=0; u<graph->nvtxs; ++u){
        random_vartex[u] = u;
    }
    srand ( time(NULL) );
    for (i = graph->nvtxs-1; i > 0; i--) {
        int j = rand() % (i+1);
        // Swap arr[i] with the element at random index
        swap_position(&random_vartex[i], &random_vartex[j]);
    }
    idx_t new_id = 0, itr = 0;
    idx_t * new_ids;
    new_ids = imalloc(graph->nvtxs, "main: part");
    int num_row = ceil(((double)graph->nvtxs)/row);
    idx_t _part = 0, nVartex = 0, start = 0;
    for (i = 0; i < graph->nvtxs; ++i) {
        new_ids[random_vartex[i]] = new_id++;
    }
    idx_t *nEdges_part = imalloc(params->nparts, "main: part");
    for (k = 0; k < params->nparts; ++k) {
        nEdges_part[k] = 0;
    }
    for (i = 0; i < graph->nvtxs; ++i) {
        nVartex++;
        for (v = graph->xadj[random_vartex[i]]; v < graph->xadj[random_vartex[i] + 1]; v++) {
            idx_t col_part = new_ids[graph->adjncy[v]]/num_row;
            nEdges_part[_part + col_part] += 1;
        }
        idx_t condition = (num_row * (_part/col + 1)) > graph->nvtxs ? graph->nvtxs : (num_row * (_part/col + 1));
        if((i+1) >= condition) {
            int startIdx = _part*num_row;
//            printf("Part=%d, row=%d, col=%d\n", _part, row, col);
            for (cl = 0; cl < col; ++cl) {
//                FILE *newMat;
//                char *ptr = strtok(params->filename, ".");
//                char outFile[MAXLINE];
//                sprintf(outFile, "%s", ptr);
//                char *ptr2 = strtok(outFile, "/");
//                char mat_filename[MAXLINE];
//                sprintf(mat_filename, "%s/output/%s_random_%"PRIDX"_%"PRIDX, ptr2, strtok(NULL, "-"), params->nparts, (_part+cl));

                FILE *newMat;
                char *last = strrchr(params->filename, '/');
                char *s = last+1;
                char *ptr = strtok(s, ".");
                char mat_filename[MAXLINE];
                sprintf(mat_filename, "random2d/%s_random_%"PRIDX"_%"PRIDX, ptr, params->nparts, (_part+cl));
                
                if (!(newMat = fopen(strcat(mat_filename, ".mtx"), "w"))) {
                    fprintf(stderr, "fopen: failed to open file '%s'", ptr);
                    exit(EXIT_FAILURE);
                }
                fprintf(newMat, "%%%MatrixMarket matrix coordinate real general\n");
                fprintf(newMat, "%d %d %d\n", nVartex, graph->nvtxs, nEdges_part[_part + cl]);

                for (itr = start; itr <= i; ++itr) {
                    u = random_vartex[itr];
                    if(new_ids[u] > i){
                        printf("%d mapped wrong as %d\n", u, new_ids[u]);
                    }
                    for (v = graph->xadj[u]; v < graph->xadj[u + 1]; v++) {
                        if(new_ids[graph->adjncy[v]]/num_row == cl) {
                            fprintf(newMat, "%d %d %lf\n", (new_ids[u] + 1), (new_ids[graph->adjncy[v]] + 1),
                                    (double) graph->adjwgt[v]);
                        }
                    }
                }
//                printf("Done=%d\n", cl);

//                 close file
                if (fclose(newMat) != 0) {
                    fprintf(stderr, "fopen: failed to open file '%s'", mat_filename);
                    exit(EXIT_FAILURE);
                }
            }
            nVartex = 0;
            start = i + 1;
            _part += col;
        }
    }
    /******* End ******/
    /***** Label the vertices with the new ID according to the partition *****/
    /*idx_t new_id = 0, itr = 0;
    idx_t * new_ids;
    idx_t * sorted_vartex;
    new_ids = imalloc(graph->nvtxs, "main: part");
    sorted_vartex = imalloc(graph->nvtxs, "main: part");
    idx_t *nVartex_part = imalloc(params->nparts, "main: part");
    idx_t *nEdges_part = imalloc(params->nparts, "main: part");
    for (k_part = 0; k_part < params->nparts; ++k_part) {
        int nVartex = 0;
        int nEdgesx = 0;
        for (u = 0; u < graph->nvtxs; u++) {
            if (part[u] == k_part) {
                sorted_vartex[new_id] = u;
                new_ids[u] = new_id++;
                nVartex++;
                if((graph->xadj[u + 1] - graph->xadj[u]) <= 0)
                    printf("No out going edges for %d\n", new_id);
                for (v = graph->xadj[u]; v < graph->xadj[u + 1]; v++) {
                    nEdgesx++;
                }
            }
        }
        nVartex_part[k_part] = nVartex;
        nEdges_part[k_part] = nEdgesx;
    }*/

    /// Convert graph into matrix in a sorting order of partition
    /*FILE *newMat;
    char *last = strrchr(params->filename, '/');
    char *s = last+1;
    char *ptr = strtok(s, ".");
    char mat_filename[MAXLINE];
    sprintf(mat_filename, "graphs/sorted/%s_%"PRIDX, ptr, params->nparts);
    if (!(newMat = fopen(strcat(mat_filename, ".mtx"), "w"))) {
        fprintf(stderr, "fopen: failed to open file '%s'", ptr);
        exit(EXIT_FAILURE);
    }
    fprintf(newMat, "%%%MatrixMarket matrix coordinate real general\n");
    fprintf(newMat, "%d %d %d\n", graph->nvtxs, graph->nvtxs, graph->nedges);
    for (k_part = 0; k_part < params->nparts; ++k_part) {
        for (itr = 0; itr < graph->nvtxs; ++itr) {
            u = sorted_vartex[itr];
            if (part[u] != k_part)
                continue;
            for (v = graph->xadj[u]; v < graph->xadj[u + 1]; v++) {
                fprintf(newMat, "%d %d %lf\n", (new_ids[u] + 1), (new_ids[graph->adjncy[v]] + 1), (double) graph->adjwgt[v]);
            }
        }
    }
    // close file
    if (fclose(newMat) != 0) {
        fprintf(stderr, "fopen: failed to open file '%s'", mat_filename);
        exit(EXIT_FAILURE);
    }*/

    /// Convert graph into matrix into multiple sorting order of partition
    int off = 1;
    if(off == 0) {
        for (k_part = 0; k_part < params->nparts; ++k_part) {
///         open file
            FILE *newMat;
            char *last = strrchr(params->filename, '/');
            char *s = last + 1;
            char *ptr = strtok(s, ".");
            char outFile[MAXLINE];
            char mat_filename[MAXLINE];
            sprintf(mat_filename, "kway/%s_%"PRIDX"_%"PRIDX, ptr, params->nparts, k_part);
            if (!(newMat = fopen(strcat(mat_filename, ".mtx"), "w"))) {
                fprintf(stderr, "fopen: failed to open file '%s'", ptr);
                exit(EXIT_FAILURE);
            }
            fprintf(newMat, "%%%MatrixMarket matrix coordinate real general\n");
            fprintf(newMat, "%d %d %d\n", nVartex_part[k_part], graph->nvtxs, nEdges_part[k_part]);

            for (itr = 0; itr < graph->nvtxs; ++itr) {
                u = sorted_vartex[itr];
                if (part[u] != k_part)
                    continue;
                for (v = graph->xadj[u]; v < graph->xadj[u + 1]; v++) {
                    fprintf(newMat, "%d %d %lf\n", (new_ids[u] + 1), (new_ids[graph->adjncy[v]] + 1),
                            (double) graph->adjwgt[v]);
                }
            }

            // close file
            if (fclose(newMat) != 0) {
                fprintf(stderr, "fopen: failed to open file '%s'", mat_filename);
                exit(EXIT_FAILURE);
            }
        }
    }
    /*FILE *nonSortMat;

    char *nonsort_ptr = strtok(params->filename, ".");
    char org_mat_filename[MAXLINE];
    sprintf(org_mat_filename, "%s", nonsort_ptr);
    if (!(nonSortMat = fopen(strcat(org_mat_filename, "_original.mtx"), "w"))) {
        fprintf(stderr, "fopen: failed to open file '%s'", nonsort_ptr);
        exit(EXIT_FAILURE);
    }
    fprintf(nonSortMat, "%%%MatrixMarket matrix coordinate real general\n");
    fprintf(nonSortMat, "%d %d %d\n", graph->nvtxs, graph->nvtxs, graph->nedges);

    for (u = 0; u < graph->nvtxs; ++u) {
        for (v = graph->xadj[u]; v < graph->xadj[u + 1]; v++) {
            fprintf(nonSortMat, "%d %d %lf\n", (u + 1), (graph->adjncy[v] + 1), (double) graph->adjwgt[v]);
        }
    }

    if (fclose(nonSortMat) != 0) {
        fprintf(stderr, "fopen: failed to open file '%s'", nonsort_ptr);
        exit(EXIT_FAILURE);
    }*/

    /// End matrix conversion

    FreeGraph(&graph);
    gk_free((void **) &part, LTERM);
    gk_free((void **) &params->filename, &params->tpwgtsfile, &params->tpwgts,
            &params->ubvecstr, &params->ubvec, &params, LTERM);

}


/*************************************************************************/
/*! This function prints run parameters */
/*************************************************************************/
void GPPrintInfo(params_t *params, graph_t *graph) {
    idx_t i;

    if (params->ufactor == -1) {
        if (params->ptype == METIS_PTYPE_KWAY)
            params->ufactor = KMETIS_DEFAULT_UFACTOR;
        else if (graph->ncon == 1)
            params->ufactor = PMETIS_DEFAULT_UFACTOR;
        else
            params->ufactor = MCPMETIS_DEFAULT_UFACTOR;
    }

    printf("******************************************************************************\n");
    printf("%s", METISTITLE);
    printf(" (HEAD: %s, Built on: %s, %s)\n", SVNINFO, __DATE__, __TIME__);
    printf(" size of idx_t: %zubits, real_t: %zubits, idx_t *: %zubits\n",
           8 * sizeof(idx_t), 8 * sizeof(real_t), 8 * sizeof(idx_t * ));
    printf("\n");
    printf("Graph Information -----------------------------------------------------------\n");
    printf(" Name: %s, #Vertices: %"PRIDX", #Edges: %"PRIDX", #Parts: %"PRIDX"\n",
           params->filename, graph->nvtxs, graph->nedges / 2, params->nparts);
    if (graph->ncon > 1)
        printf(" Balancing constraints: %"PRIDX"\n", graph->ncon);

    printf("\n");
    printf("Options ---------------------------------------------------------------------\n");
    printf(" ptype=%s, objtype=%s, ctype=%s, rtype=%s, iptype=%s\n",
           ptypenames[params->ptype], objtypenames[params->objtype], ctypenames[params->ctype],
           rtypenames[params->rtype], iptypenames[params->iptype]);

    printf(" dbglvl=%"PRIDX", ufactor=%.3f, no2hop=%s, minconn=%s, contig=%s, nooutput=%s\n",
           params->dbglvl,
           I2RUBFACTOR(params->ufactor),
           (params->no2hop ? "YES" : "NO"),
           (params->minconn ? "YES" : "NO"),
           (params->contig ? "YES" : "NO"),
           (params->nooutput ? "YES" : "NO")
    );

    printf(" seed=%"PRIDX", niter=%"PRIDX", ncuts=%"PRIDX"\n",
           params->seed, params->niter, params->ncuts);

    if (params->ubvec) {
        printf(" ubvec=(");
        for (i = 0; i < graph->ncon; i++)
            printf("%s%.2e", (i == 0 ? "" : " "), (double) params->ubvec[i]);
        printf(")\n");
    }

    printf("\n");
    switch (params->ptype) {
        case METIS_PTYPE_RB:
            printf("Recursive Partitioning ------------------------------------------------------\n");
            break;
        case METIS_PTYPE_KWAY:
            printf("Direct k-way Partitioning ---------------------------------------------------\n");
            break;
    }
}


/*************************************************************************/
/*! This function does any post-partitioning reporting */
/*************************************************************************/
void GPReportResults(params_t *params, graph_t *graph, idx_t *part, idx_t objval) {
    gk_startcputimer(params->reporttimer);
    ComputePartitionInfo(params, graph, part);

    gk_stopcputimer(params->reporttimer);

    printf("\nTiming Information ----------------------------------------------------------\n");
    printf("  I/O:          \t\t %7.3"PRREAL" sec\n", gk_getcputimer(params->iotimer));
    printf("  Partitioning: \t\t %7.3"PRREAL" sec   (METIS time)\n", gk_getcputimer(params->parttimer));
    printf("  Reporting:    \t\t %7.3"PRREAL" sec\n", gk_getcputimer(params->reporttimer));
    printf("\nMemory Information ----------------------------------------------------------\n");
    printf("  Max memory used:\t\t %7.3"PRREAL" MB\n", (real_t)(params->maxmemory / (1024.0 * 1024.0)));
    printf("******************************************************************************\n");

}
