#ifndef LSDEBUG_H
#define LSDEBUG_H
#define MYDBG printf("%s:%d\n",__FILE__,__LINE__);
#define MYDBG0 if(cf.rank==0) printf("%s:%d\n",__FILE__,__LINE__);
#endif
