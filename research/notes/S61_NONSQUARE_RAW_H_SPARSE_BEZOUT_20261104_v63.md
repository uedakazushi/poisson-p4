# Raw H14×5の固定四行・一列だけのQQ係数span

固定native行2,4,5,8は残りtargetの0と5にだけsupportを持つ。一列0の四係数
だけを使い、既知a,d,F,B,C,J因子と非零QQ定数を除いたprimitiveに対し、定数、
次にaffine(1,a,e,d) cofactorを試す。未知数は4又は16。targetは厳密に
1,a,d,F,B,C,J,adの八つだけで、QQ係数の有限連立一次方程式を解く。
点での評価、GB、一般のminor探索は行わない。

成功した場合はprimitive恒等式だけで終えず、元四行の係数へ既知単元を戻し、
実際の4×2行の全成分を確認する。他方のcolumn5はtailとしてそのまま残す。
さらに前段のrow清掃と一つのpivot、Ri=Ci−alpha_i*Cdeltaを逆に辿り、元genuine16
のfull source係数を保存する。全7個の実際のH係数を元16×7から再構成して検査
する。pure cubic部分は同じ16係数を使うlazy式として保持する。

単一columnの単元と非零tailを得てもsource基底のcompletionやtarget消去を
行ったことにはしない。成功・不成功のいずれも限定したcofactor/target classの
結果として保存し、fullH7/full57/full239の結論へ飛躍しない。
