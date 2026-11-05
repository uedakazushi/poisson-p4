# Dense known-adjoint worker accounting

The proposed worker stores normalized210-coordinate vectors and their91-coordinate original-source representations as dense coefficient arrays. Original sparse source columns remain shared input references. This avoids per-entry numeric index duplication and includes every pivot-step record, two current vector copies, target representations, phase metadata, scalar scratch and one-row serialization overlap.

The explicit whole bound is51562 numeric slots within the unchanged58603 cap. A single row-streamed atomic frontier, its temporary predecessor, one append-only pivot history, input bytes and result metadata use at most16842088 bytes within the unchanged16940016 cap. Full frontiers are not accumulated. The old sparse-method expression is not a whole-worker proof; no solve ran under it. This distinct clarification changes no mathematical input or resource cap.

The384-bit committed and1153-bit uncommitted limits are dynamic stopping guards, not an a priori completed-rank-height theorem. A failure must preserve the last complete basis, representations, current residual, pivot history and pending operation without retry. There is no physical, gauge or unknown coefficient evaluation.
