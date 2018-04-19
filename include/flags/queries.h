
/*******************************************************************************
*
* File flags/queries.h
*
* Copyright (C) 2009-2012, 2016 Martin Luescher, Isabel Campos
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Query descriptions
*
*******************************************************************************/

#define QUERIES_H

#if (defined FLAGS_C)

static int (*query_fcts[(int)(QUERIES)+1])(void)={NULL};

static int QueryUMatchUd(void)
{
   return (lat.u==lat.ud);
}

static int QueryUdbufUp2date(void)
{
   return ((lat.ud>0)&&(lat.udbuf==lat.ud));
}

static int QueryBstapUp2date(void)
{
   return ((lat.ud>0)&&(lat.bstap==lat.ud));
}

static int QueryFtsUp2date(void)
{
   return ((lat.ud>0)&&(lat.fts==lat.ud));
}

static int QuerySwUp2date(void)
{
   return ((lat.u>0)&&(lat.sw[0]==lat.u));
}

static int QuerySwEInverted(void)
{
   return (lat.sw[1]==1);
}

static int QuerySwOInverted(void)
{
   return (lat.sw[2]==1);
}

static int QuerySwdUp2date(void)
{
   return ((lat.ud>0)&&(lat.swd[0]==lat.ud));
}

static int QuerySwdEInverted(void)
{
   return (lat.swd[1]==1);
}

static int QuerySwdOInverted(void)
{
   return (lat.swd[2]==1);
}

static int QueryAwUp2date(void)
{
   return ((lat.ud>0)&&(lat.aw==lat.ud));
}

static int QueryAwhatUp2date(void)
{
   return ((lat.ud>0)&&(lat.awh==lat.ud));
}

static int QueryUdPhaseSet(void)
{
   return (lat.phase==1);
}

static void set_queries(void)
{
   query_fcts[(int)(U_MATCH_UD)]=QueryUMatchUd;
   query_fcts[(int)(UDBUF_UP2DATE)]=QueryUdbufUp2date;
   query_fcts[(int)(BSTAP_UP2DATE)]=QueryBstapUp2date;
   query_fcts[(int)(FTS_UP2DATE)]=QueryFtsUp2date;
   query_fcts[(int)(SW_UP2DATE)]=QuerySwUp2date;
   query_fcts[(int)(SW_E_INVERTED)]=QuerySwEInverted;
   query_fcts[(int)(SW_O_INVERTED)]=QuerySwOInverted;
   query_fcts[(int)(SWD_UP2DATE)]=QuerySwdUp2date;
   query_fcts[(int)(SWD_E_INVERTED)]=QuerySwdEInverted;
   query_fcts[(int)(SWD_O_INVERTED)]=QuerySwdOInverted;
   query_fcts[(int)(AW_UP2DATE)]=QueryAwUp2date;
   query_fcts[(int)(AWHAT_UP2DATE)]=QueryAwhatUp2date;
   query_fcts[(int)(UD_PHASE_SET)]=QueryUdPhaseSet;
}

#endif
