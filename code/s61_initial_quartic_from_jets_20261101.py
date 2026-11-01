"""Exact block construction of I_2,I_3,I_4 from one consistent raw 4-jet.
Input rows are homogeneous coefficient matrices A2,A3,A4 in lexicographic
combinations_with_replacement monomial order. Coefficient domain is a FIELD.
No assertion about degrees above4 or the full initial ideal is made.
"""
from itertools import combinations_with_replacement
import sympy as s
from sympy.polys.matrices import DomainMatrix


def initial_quartic(A2,A3,A4,n,domain,compute_rank=False,progress=lambda x:None):
    assert domain.is_Field, "The coefficient domain must be a field."
    matrices=[s.Matrix(A) for A in [A2,A3,A4]]
    count=matrices[0].rows
    mons={d:list(combinations_with_replacement(range(n),d))for d in [2,3,4]}
    assert all(A.rows==count and A.cols==len(mons[d])for d,A in zip([2,3,4],matrices))
    a2,a3,a4=[A.to_DM().convert_to(domain)for A in matrices]
    def extract(A,rows,cols):return A.extract(list(rows),list(cols))
    def independent_rows(A):return list(A.transpose().rref()[1])
    def express(A,B):
        piv=list(B.rref()[1]);assert len(piv)==B.shape[0]
        C=extract(A,range(A.shape[0]),piv)*extract(B,range(B.shape[0]),piv).inv()
        assert C*B==A
        return C,piv
    def multiples(A,degree):
        source=A.to_Matrix();idx={m:i for i,m in enumerate(mons[degree+1])};out={}
        for row in range(source.rows):
            for j in range(n):
                for k,m in enumerate(mons[degree]):
                    c=source[row,k]
                    if c:out[row*n+j,idx[tuple(sorted(m+(j,)))]]=c
        return s.SparseMatrix(source.rows*n,len(mons[degree+1]),out).to_DM().convert_to(domain)
    rows2=independent_rows(a2);b2=extract(a2,rows2,range(a2.shape[1]));c2,piv2=express(a2,b2)
    f3=extract(a3,rows2,range(a3.shape[1]));f4=extract(a4,rows2,range(a4.shape[1]))
    h3=a3-c2*f3;h4=a4-c2*f4
    lead3=multiples(b2,2).vstack(h3)
    tail4=multiples(f3,3).vstack(h4)
    rows3=independent_rows(lead3);b3=extract(lead3,rows3,range(lead3.shape[1]));c3,piv3=express(lead3,b3)
    progress(dict(quadratic_rank=len(rows2),cubic_initial_rank=len(rows3),degree3_source_rows=lead3.shape[0]))
    extra4=tail4-c3*extract(tail4,rows3,range(tail4.shape[1]))
    lower4=multiples(b3,3)
    result=dict(quadratic_rank=len(rows2),cubic_initial_rank=len(rows3),quadratic_basis_rows=rows2,cubic_basis_source_rows=rows3,quadratic_pivots=piv2,cubic_pivots=piv3,quadratic_row_reduction=c2,cubic_row_reduction=c3,degree3_leading=lead3,degree3_tail=tail4,quartic_lower=lower4,quartic_extra=extra4)
    if compute_rank:
        result['quartic_lower_rank']=len(lower4.rref()[1])
        result['quartic_initial_rank']=len(lower4.vstack(extra4).rref()[1])
    return result
