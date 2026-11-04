# G115の第十二pivot：gとD5の追加chartを残さない3対10

2026-11-04。既計算N6の係数を再展開せず、source-only有限operatorを作る。
N6は固定六minorの正規化値で9047項非零、元D6=a二乗N6である。
旧四対十一のG115/target19成分をrとすると、既存比例関係から

    g D5_G20 r = c10 a² N6

が成り立つ。初めはD(gD5)でborder比から得たが、両辺は旧
FULLU[Delta407逆,Saux逆,T逆]の元であり、整域の稠密開集合上の一致から
この局所環そのもので一致する。

従ってN6を単元にした開集合ではg、D5、rも全て単元であり、

    r逆 = g D5_G20 /(c10 a² N6)

を使える。g逆やD5逆という新chartは不要である。D5はliteral source五行と
target五列の有限行列式recipeで保持し、係数は計算しない。全補集合、特にN6零集合は保持する。

旧residual sourceの第二行G115と最後のtarget19を先頭へ動かすと
[r B;C D]となる。左[rinv 0;−C rinv I]と右[1 −rinv B;0 I]を掛け、
diag(1,D−C rinv B)とする。両側の明示逆を旧15source/22target逆へ合成し、
FULLU D(Delta407 Saux T N6)上の同cokernel三対十を得る。

残るsourceはcurrent15の[9,12,14]、targetはoriginal32の
[1,3,4,6,7,8,9,14,17,18]である。元八追加tailのうち既消費2以外の七tail効果を
全て保持する。三sourceのold25、normalized full57、canonical63への係数を有限operatorで
明示し、各sourceのHZ/HWをeta単元とCdeltaで修正して実raw七heavy列を零にする。
ここでraw未正規化Jacobi135との係数同一性を新たに主張しない。

この三対十がrank三なら元cubic57が成立する、という十分条件を保持する。
三行のrank、N6零集合上の階数、全FULLU57とquartic239および尽くしは未決である。
