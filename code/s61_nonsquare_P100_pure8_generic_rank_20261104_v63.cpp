// Incremental selected square only; never eliminate the full8x98 matrix.
#define main frozen_q5_row_main_NOT_CALLED
#include "s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp"
#undef main
bool iszero(const P&p){return fmpq_poly_is_zero(p.v);}
bool same(const P&a,const P&b){return fmpq_poly_equal(a.v,b.v);}
P minuspoly(const P&a,const P&b){return plusp(a,negp(b));}
using Mat=vector<vector<P>>;
void checkadj(const Mat&S,const Mat&H,const P&D){
 int n=S.size();assert(!iszero(D));for(int i=0;i<n;i++)for(int j=0;j<n;j++){P l,r;for(int k=0;k<n;k++){l=plusp(l,mulp(S[i][k],H[k][j]));r=plusp(r,mulp(H[i][k],S[k][j]));}P want=i==j?D:P();assert(same(l,want)&&same(r,want));}
}
void emitmatrix(ostream&o,const Mat&A){o<<'[';for(int i=0;i<int(A.size());i++){if(i)o<<',';o<<'[';for(int j=0;j<int(A[i].size());j++){if(j)o<<',';jp(o,A[i][j]);}o<<']';}o<<']';}
void emitints(ostream&o,const vector<int>&v){o<<'[';for(int i=0;i<int(v.size());i++){if(i)o<<',';o<<v[i];}o<<']';}
int main(int argc,char**argv){assert(argc==3);began=chrono::steady_clock::now();string folder=argv[2];try{
 ifstream f(argv[1]);f>>input_reserve;vector<vector<int>>ex(8,vector<int>(5));for(auto&row:ex)for(auto&v:row)f>>v;
 P zz=shiftp(one(),1);vector<P>factors={zz,plusp(zz,pc("20")),plusp(zz,pc("60")),plusp(scalar(zz,"203"),pc("2700")),plusp(zz,pc("12"))};
 Mat A(8,vector<P>(98));vector<P>LCD(8);stage="all784_known_unit_rowLCD_clearings";
 for(int i=0;i<8;i++){P den=one();for(int k=0;k<5;k++)for(int j=0;j<ex[i][k];j++)den=mulp(den,factors[k]);LCD[i]=den;
  for(int j=0;j<98;j++){R r=inputr(f);P q=quotient(den,r.d);A[i][j]=mulp(r.n,q);assert(same(mulp(A[i][j],r.d),mulp(r.n,den)));}
 }string extra;assert(!(f>>extra));
 Mat S,H;P D=one();vector<int>rows,cols;vector<P>detsteps;long candidate_count=0,division_count=0;stage="incremental_exact_selected_minor";
 for(int column=0;column<98&&rows.size()<8;column++){
  int n=rows.size();vector<P>b(n),u(n);for(int j=0;j<n;j++)b[j]=A[rows[j]][column];for(int i=0;i<n;i++)for(int j=0;j<n;j++)u[i]=plusp(u[i],mulp(H[i][j],b[j]));
  for(int row=0;row<8;row++){
   if(find(rows.begin(),rows.end(),row)!=rows.end())continue;candidate_count++;vector<P>c(n),v(n);for(int j=0;j<n;j++)c[j]=A[row][cols[j]];
   P h=mulp(D,A[row][column]);for(int j=0;j<n;j++)h=minuspoly(h,mulp(c[j],u[j]));if(iszero(h))continue;
   for(int j=0;j<n;j++)for(int k=0;k<n;k++)v[j]=plusp(v[j],mulp(c[k],H[k][j]));
   Mat NS(n+1,vector<P>(n+1)),NH(n+1,vector<P>(n+1));
   for(int i=0;i<n;i++)for(int j=0;j<n;j++){NS[i][j]=S[i][j];NH[i][j]=quotient(plusp(mulp(h,H[i][j]),mulp(u[i],v[j])),D);division_count++;}
   for(int i=0;i<n;i++){NS[i][n]=b[i];NS[n][i]=c[i];NH[i][n]=negp(u[i]);NH[n][i]=negp(v[i]);}
   NS[n][n]=A[row][column];NH[n][n]=D;checkadj(NS,NH,h);rows.push_back(row);cols.push_back(column);S=move(NS);H=move(NH);D=h;detsteps.push_back(D);break;
  }
 }
 int rank=rows.size();assert(rank>0&&!iszero(D));Mat C;
 if(rank<8){stage="all784_row_span_upper_identities";C=Mat(8,vector<P>(rank));for(int i=0;i<8;i++)for(int j=0;j<rank;j++)for(int k=0;k<rank;k++)C[i][j]=plusp(C[i][j],mulp(A[i][cols[k]],H[k][j]));
  for(int i=0;i<8;i++)for(int j=0;j<98;j++){P v;for(int k=0;k<rank;k++)v=plusp(v,mulp(C[i][k],A[rows[k]][j]));assert(same(v,mulp(D,A[i][j])));}
 }
 checkadj(S,H,D);stage="compact_witness_serialization";
 {ofstream o(folder+"/native.json");o<<"{\"generic_rank\":"<<rank<<",\"all784_known_unit_row_clearings_verified\":true,\"selected_rows\":";emitints(o,rows);o<<",\"selected_columns\":";emitints(o,cols);o<<",\"determinant_nonzero\":";jp(o,D);o<<",\"selected_square\":";emitmatrix(o,S);o<<",\"final_adjugate\":";emitmatrix(o,H);o<<",\"all_selected_adjugate_left_right_products_verified\":true,\"exact_division_identity_count\":"<<division_count<<",\"determinant_steps\":[";for(int i=0;i<rank;i++){if(i)o<<',';jp(o,detsteps[i]);}o<<"],\"row_span_coordinates_numerator8xr\":";emitmatrix(o,C);o<<",\"all784_span_upper_identities_verified\":"<<(rank<8?"true":"false")<<",\"rank_upper_from_number_of_rows_if8\":"<<(rank==8?"true":"false")<<",\"candidate_bordered_entries_checked\":"<<candidate_count<<",\"peak_live_polynomial_terms_with_literal_input\":"<<peak_live+input_reserve<<",\"serialized_terms\":"<<stored<<",\"peak_polynomial_length\":"<<peak<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";}
 }catch(const exception&e){ofstream o(folder+"/frontier.json");o<<"{\"all_passed\":false,\"stage\":\""<<stage<<"\",\"reason\":\""<<e.what()<<"\"}\n";return 2;}flint_cleanup();}
