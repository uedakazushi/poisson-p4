# 固定K120/target14分子のnested一回計算

2026-11-04。独立source gate 2c48439bの後、同じ分子について一回だけ実行する。
原表示の三項を保持し、gは数値除法せず形式交代性で先に除く。
source別のcurrent15行を用い、bare Jacobi135行とは同一視しない。

各外側項でA/B/Cの二列だけのweighted二minorを計算し、それ以外のsource cacheを解放する。
次に内側整数格子を確保し、十個の二minor×三minorをそれぞれ新しいsource cacheから
作って融合する。唯一の外側整数格子へ融合した後、内側とweighted値を全て解放する。
全210内scale・21外scaleを理論分母と実contentの両方で検査する。

global atom、phase source値、二整数格子、scratchと15000作業枠を同時に数え、
150kを越えたら止める。内ループの時計guardは五十秒で、部分停止なら完了した
outer/inner数、現在source節点と畳み込み位置、正確な二整数格子を保存する。
停止prefixを全分子や非零判定と呼ばず、同じworkerの自動再試行もしない。
完了時も一つのglobal分母と整数係数列で完全な有理多項式を保存する。

生のD7/gは既知a二乗×出力であり、現在の残余一成分はc11×出力/N6となる。
この一成分の結果を他の列・残余rank・全FULLU57へ拡張しない。
全tail2効果、補集合、以前のraw十minorと通常格子の費用失敗を保持する。
