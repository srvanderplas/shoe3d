#' Read in a shoe from a file path
#'
#' @param filepath path to shoe STL file
#' @return centered shoe mesh object
#' @export
#' @importFrom Rvcg vcgImport
#' @importFrom Morpho barycenter
#' @importFrom rgl translate3d
#' @import assertthat
input_shoe <- function(filepath) {

  # checking to see if the file exists
  assertthat::see_if(file.exists(filepath))

  # turning it into a mesh object
  shoe_mesh <- Rvcg::vcgImport(filepath, clean = T)
  # Checking that it is a mesh3d object
  assertthat::assert_that(class(shoe_mesh) == "mesh3d")
  centering <- Morpho::barycenter(shoe_mesh) %>% colMeans()
  shoemesh <- rgl::translate3d(shoe_mesh, -centering[1], -centering[2], -centering[3])

  return(shoemesh)
}

#' Get coordinates of mesh vertices
#'
#' This fuction turns a mesh object into a matrix of coordinates
#'
#' @param shoe the mesh object of a shoe from the data set
#' @param verts the number of edges to each vertex (if na does full shoe) to filter out of the shoe
#' @return data.frame
#'
#'
#' @importFrom geometry convhulln
#' @importFrom Arothron trasf.mesh bary.mesh
#' @importFrom tidyr gather spread
#' @import dplyr
#' @import assertthat

shoe_coord <- function(shoe, verts = NULL) {
  if (!is.null(verts)) {
    vert <- as.data.frame(t(shoe$it)) %>%
      dplyr::mutate(triangle_id = 1:n()) %>%
      tidyr::gather(-triangle_id, key = "vertex", value = idx) %>%
      dplyr::group_by(idx) %>%
      dplyr::filter(n() <= verts) %>%
      dplyr::ungroup(idx) %>%
      tidyr::spread(key = "vertex", value = idx) %>%
      dplyr::filter(!is.na(V1) & !is.na(V2) & !is.na(V3))
  } else {
    vert <- as.data.frame(t(shoe$it)) %>%
      dplyr::mutate(triangle_id = 1:n()) %>%
      tidyr::gather(-triangle_id, key = "vertex", value = idx) %>%
      tidyr::spread(key = "vertex", value = idx) %>%
      dplyr::filter(!is.na(V1) & !is.na(V2) & !is.na(V3))
  }

  vert_long <- vert %>%
    dplyr::select(-triangle_id) %>%
    tidyr::gather(key = V, value = idx) %>%
    dplyr::select(-V) %>%
    unique()
  vert_coords <- shoe %>%
    `[[`("vb") %>%
    t() %>%
    magrittr::set_colnames(c("x", "y", "z", "idk")) %>%
    as.data.frame() %>%
    dplyr::mutate(idx = 1:n())


  x <- dplyr::left_join(vert_long, vert_coords) %>%
    dplyr::select(-idk) %>%
    dplyr::select(-idx)

  hul <- geometry::convhulln(x,
    output.options = c("p", "Fx"),
    return.non.triangulated.facets = FALSE
  )



  hull <- as.vector(hul) %>%
    unique() %>%
    as.data.frame()
  colnames(hull) <- "idx"

  edge_coords <- vert_coords %>%
    dplyr::filter(idx %in% hull$idx) %>%
    select(-idk) %>%
    select(-idx) %>%
    as.matrix()

  # assertthat::assert_that(class(edge_coords) == "matrix")

  return(edge_coords)
}

#' Do a simple check to see if rotations match by minimizing RMSE
#'
#' @param m1 mesh 1
#' @param m2 mesh 2
#' @return mesh1 oriented in the same position as mesh 2.
#'
#' @export

fix_rotation <- function(m1, m2) {
  rotations <- list(
    rotationMatrix(matrix=diag(c(1,1,1))),
    rotationMatrix(matrix=diag(c(-1,1,1))),
    rotationMatrix(matrix=diag(c(1,-1,1))),
    rotationMatrix(matrix=diag(c(1,1,-1))),
    rotationMatrix(matrix=diag(c(1,-1,-1))),
    rotationMatrix(matrix=diag(c(-1,1,-1))),
    rotationMatrix(matrix=diag(c(-1,-1,1))),
    rotationMatrix(matrix=diag(c(-1,-1,-1)))
  )

  m1_trans <- purrr::map(rotations, ~transform3d(m1, .))
  dists <- purrr::map(m1_trans, ~Morpho::meshDist(., m2, plot = F))
  rmse <- purrr::map_dbl(dists, ~sqrt(mean(.$dists^2)))

  return(m1_trans[which.min(rmse)])
}
