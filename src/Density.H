namespace transport_methods {
enum transport_method {
  TEST=0,
  NEGF=1,
  WF=2,
};
} // namespace

int density(TCSR<double> *,TCSR<double> *,TCSR<double> *,CPX,CPX,transport_methods::transport_method,c_transport_type);
