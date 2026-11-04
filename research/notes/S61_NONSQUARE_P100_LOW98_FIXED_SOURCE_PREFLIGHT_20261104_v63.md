# P100 全239への fixed low-source frontier

既受理 H15 の七剰余体の full239 証明を source-only で比較した。七つの証明は同じ順序の98個の original319 source 行を選ぶ。その内訳は quadratic product11行と cubic product87行で、genuine quartic source は含まない。quadratic 行は `[0,1,2,4,5,7,8,10,11,12,145]`。元354→319の保存ラベルで照合した。

従って七つの別々の source 行列を持ち上げる必要はなく、一つの actual global lazy matrix B(z) と determinant D(z) を指定できる。有限体係数を任意に有理数へ持ち上げる意味ではない。D の全七因子上での非零性は既受理証明から従う。D は非零な実 minor であり、D の非零開では full239 が成立する。

ただし D の全係数・因数・単元性は未計算である。H15 上で D が単元でも、その外の有限 torsion cokernel は除外できない。pure89 の同一 cokernel 表示も、この有限 torsion を自動的には消さない。曲線の PID 上で generic kernel が自由であることだけでは、fiberwise kernel の次元増加を防げない。

この比較自体は新しい scalar 係数、normal form、rank、minor、点を計算しない。全域239は未決である。次の実計算が必要なら、この固定 low-source witness を source anchor として用い、実 global coefficient/support bound と単元 ideal の証明を別に行う。98×98を無条件に展開する指示ではない。
