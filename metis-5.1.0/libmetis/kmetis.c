/*!
\file  
\brief The top-level routines for  multilevel k-way partitioning that minimizes
       the edge cut.

\date   Started 7/28/1997
\author George  
\author Copyright 1997-2011, Regents of the University of Minnesota 
\version\verbatim $Id: kmetis.c 13905 2013-03-25 13:21:20Z karypis $ \endverbatim
*/

#include "metislib.h"


/*************************************************************************/
/*! This function is the entry point for MCKMETIS */
/*************************************************************************/
int METIS_PartGraphKway(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy, 
          idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *nparts, 
          real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *objval, 
          idx_t *part)
{
//    printf("Start k way partition\n");
  int sigrval=0, renumber=0;
  graph_t *graph;
  ctrl_t *ctrl;

  /* set up malloc cleaning code and signal catchers */
  if (!gk_malloc_init()) 
    return METIS_ERROR_MEMORY;
  gk_sigtrap();

  if ((sigrval = gk_sigcatch()) != 0)
    goto SIGTHROW;

  /* set up the run parameters */
  ctrl = SetupCtrl(METIS_OP_KMETIS, options, *ncon, *nparts, tpwgts, ubvec);
  if (!ctrl) {
      printf("setup ctrl failed\n");
    gk_siguntrap();
    return METIS_ERROR_INPUT;
  }
  /* if required, change the numbering to 0 */
  if (ctrl->numflag == 1) {
    Change2CNumbering(*nvtxs, xadj, adjncy);
    renumber = 1;
  }
  /* set up the graph */
  if(vsize == NULL){
      vsize = ismalloc((*nvtxs), 1, "ReadGraph: vsize");
  }
  if(vwgt == NULL){
      printf("vwgt null found\n");
      vwgt = ismalloc((*ncon)*(*nvtxs), 1, "ReadGraph: vwgt");
  }
    graph = SetupGraph(ctrl, *nvtxs, *ncon, xadj, adjncy, vwgt, vsize, adjwgt);
    /* set up multipliers for making balance computations easier */
  SetupKWayBalMultipliers(ctrl, graph);
  /* set various run parameters that depend on the graph */
  ctrl->CoarsenTo = gk_max((*nvtxs)/(20*gk_log2(*nparts)), 30*(*nparts));
  printf("ctrl->CoarsenTo %d from %d and %d\n", ctrl->CoarsenTo, ((*nvtxs)/(20*gk_log2(*nparts))), (30*(*nparts)));
  ctrl->nIparts   = (ctrl->CoarsenTo == 30*(*nparts) ? 4 : 5);
  printf("ctrl->nIparts %d\n", ctrl->nIparts);
  /* take care contiguity requests for disconnected graphs */
  if (ctrl->contig && !IsConnected(graph, 0)) 
    gk_errexit(SIGERR, "METIS Error: A contiguous partition is requested for a non-contiguous input graph.\n");
  /* allocate workspace memory */  
  AllocateWorkSpace(ctrl, graph);
  /* start the partitioning */
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->TotalTmr));
  *objval = MlevelKWayPartitioning(ctrl, graph, part);
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->TotalTmr));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));
  /* clean up */
  FreeCtrl(&ctrl);
SIGTHROW:
  /* if required, change the numbering back to 1 */
  if (renumber)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
  gk_siguntrap();
    gk_malloc_cleanup(0);
    return metis_rcode(sigrval);
}


/*************************************************************************/
/*! This function computes a k-way partitioning of a graph that minimizes
    the specified objective function.

    \param ctrl is the control structure
    \param graph is the graph to be partitioned
    \param part is the vector that on return will store the partitioning

    \returns the objective value of the partitoning. The partitioning 
             itself is stored in the part vector.
*/
/*************************************************************************/
idx_t MlevelKWayPartitioning(ctrl_t *ctrl, graph_t *graph, idx_t *part)
{
  idx_t i, j, objval=0, curobj=0, bestobj=0;
  real_t curbal=0.0, bestbal=0.0;
  graph_t *cgraph;
  int status;

//printf("ctrl->ncuts %d\n", ctrl->ncuts);
  for (i=0; i<ctrl->ncuts; i++) {
//      printf("0\n");
    cgraph = CoarsenGraph(ctrl, graph);
//      printf("1\n");
    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->InitPartTmr));
    AllocateKWayPartitionMemory(ctrl, cgraph);
//printf("2\n");
    /* Release the work space */
    FreeWorkSpace(ctrl);
//      printf("3\n");
    /* Compute the initial partitioning */
    InitKWayPartitioning(ctrl, cgraph);
//      printf("4\n");
    /* Re-allocate the work space */
    AllocateWorkSpace(ctrl, graph);
    AllocateRefinementWorkSpace(ctrl, 2*cgraph->nedges);
//      printf("5\n");
    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->InitPartTmr));
    IFSET(ctrl->dbglvl, METIS_DBG_IPART, 
        printf("Initial %"PRIDX"-way partitioning cut: %"PRIDX"\n", ctrl->nparts, objval));
//      printf("6\n");
    RefineKWay(ctrl, graph, cgraph);
//      printf("7\n");
    switch (ctrl->objtype) {
      case METIS_OBJTYPE_CUT:
        curobj = graph->mincut;
        break;

      case METIS_OBJTYPE_VOL:
        curobj = graph->minvol;
        break;

      default:
        gk_errexit(SIGERR, "Unknown objtype: %d\n", ctrl->objtype);
    }
//      printf("8\n");
    curbal = ComputeLoadImbalanceDiff(graph, ctrl->nparts, ctrl->pijbm, ctrl->ubfactors);
//      printf("9\n");
    if (i == 0 
        || (curbal <= 0.0005 && bestobj > curobj)
        || (bestbal > 0.0005 && curbal < bestbal)) {
      icopy(graph->nvtxs, graph->where, part);
      bestobj = curobj;
      bestbal = curbal;
    }
//      printf("10\n");
    FreeRData(graph);
//      printf("11\n");
    if (bestobj == 0)
      break;
  }
//    printf("12\n");
  FreeGraph(&graph);
//    printf("13\n");
  return bestobj;
}


/*************************************************************************/
/*! This function computes the initial k-way partitioning using PMETIS 
*/
/*************************************************************************/
void InitKWayPartitioning(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, ntrials, options[METIS_NOPTIONS], curobj=0, bestobj=0;
  idx_t *bestwhere=NULL;
  real_t *ubvec=NULL;
  int status;

  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NITER]   = 10;
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_NO2HOP]  = ctrl->no2hop;


  ubvec = rmalloc(graph->ncon, "InitKWayPartitioning: ubvec");
  for (i=0; i<graph->ncon; i++) 
    ubvec[i] = (real_t)pow(ctrl->ubfactors[i], 1.0/log(ctrl->nparts));


  switch (ctrl->objtype) {
    case METIS_OBJTYPE_CUT:
    case METIS_OBJTYPE_VOL:
      options[METIS_OPTION_NCUTS] = ctrl->nIparts;
      status = METIS_PartGraphRecursive(&graph->nvtxs, &graph->ncon, 
                   graph->xadj, graph->adjncy, graph->vwgt, graph->vsize, 
                   graph->adjwgt, &ctrl->nparts, ctrl->tpwgts, ubvec, 
                   options, &curobj, graph->where);

      printf("METIS_PartGraphRecursive\n");
      if (status != METIS_OK)
        gk_errexit(SIGERR, "Failed during initial partitioning\n");

      break;

#ifdef XXX /* This does not seem to help */
      print("test does this code execute or not \n");
    case METIS_OBJTYPE_VOL:
      bestwhere = imalloc(graph->nvtxs, "InitKWayPartitioning: bestwhere");
      options[METIS_OPTION_NCUTS] = 2;

      ntrials = (ctrl->nIparts+1)/2;
      for (i=0; i<ntrials; i++) {
        status = METIS_PartGraphRecursive(&graph->nvtxs, &graph->ncon, 
                     graph->xadj, graph->adjncy, graph->vwgt, graph->vsize, 
                     graph->adjwgt, &ctrl->nparts, ctrl->tpwgts, ubvec, 
                     options, &curobj, graph->where);
        if (status != METIS_OK)
          gk_errexit(SIGERR, "Failed during initial partitioning\n");

        curobj = ComputeVolume(graph, graph->where);

        if (i == 0 || bestobj > curobj) {
          bestobj = curobj;
          if (i < ntrials-1)
            icopy(graph->nvtxs, graph->where, bestwhere);
        }

        if (bestobj == 0)
          break;
      }
      if (bestobj != curobj)
        icopy(graph->nvtxs, bestwhere, graph->where);

      break;
#endif

    default:
      gk_errexit(SIGERR, "Unknown objtype: %d\n", ctrl->objtype);
  }

  gk_free((void **)&ubvec, &bestwhere, LTERM);

}


