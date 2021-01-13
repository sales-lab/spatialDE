.set_fake_tqdm <- "
def fake_tqdm(iterable=None, *args, **kwargs):
  assert iterable is not None
  yield from iterable

import SpatialDE
SpatialDE.base.tqdm = fake_tqdm
"

.set_real_tqdm <- "
import SpatialDE
import tqdm.autonotebook
SpatialDE.base.tqdm = tqdm.autonotebook.tqdm
"

#' Import SpatialDE
#'
#' This function loads the SpatialDE Python module and optionally
#' monkey-patches it to remove tqdm calls.
#'
#' @param patch_tqdm If `TRUE` patch calls to tqdm.
#' @return An R wrapper for the SpatialDE Python module.
#'
#' @importFrom reticulate import py_run_string
#'
.importPyModule <- function(patch_tqdm) {
  mod <- import("SpatialDE")
  py_run_string(ifelse(patch_tqdm, .set_fake_tqdm, .set_real_tqdm))
  return(mod)
}
