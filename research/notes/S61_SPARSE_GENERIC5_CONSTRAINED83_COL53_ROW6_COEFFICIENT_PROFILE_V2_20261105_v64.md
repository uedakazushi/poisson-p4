# 固定row6/col53 SOURCE profile V2（prepared）

実行仕様は同日COL53_ROW6_PROFILE_V2_SOURCE_GATEを参照。旧失敗worker/helperはそのまま残し、新しいhelper/controllerにだけ一走査card/BOX/Hcounts・境界hash・全LIVE snapshotを実装する。SOURCE gate SHAとROOTレビューSHAの両方を固定した --profile 以外は実支持計算を始めない。SOURCEの関数所有者・式・raw11成分・旧8境界は変えない。

profile未実行時点のprepared noteであり、全11成分や全14行の完成を主張しない。係数VALUE workerのGOでもない。将来のH区間によるVALUE backpropには、その全途中支持・codec・storage・readbackを別に証明する必要がある。
