// contains the parts common to the sbrueck/mcalderara codes later to be merged

#ifdef __cplusplus
extern "C" {
#endif

   typedef void(*ptr2function) (c_DBCSR, c_DBCSR, c_DBCSR *);

   void c_scf_method(c_DBCSR S, c_DBCSR KS, c_DBCSR * P);

   void yyrestart(FILE *);
   void yyparse();
   extern FILE *yyin;

#ifdef __cplusplus
}
#endif
