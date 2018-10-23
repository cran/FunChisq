# plot-table.R
#
# Joe Song
# Created: October adapted from discrete-patterns.R April 3, 2017

plot_table <- function(table, xlab="Column", ylab="Row", col="green3",
                       xaxt="n", yaxt="n", main=NA, show.value=TRUE, ...)
{
  if(!is.matrix(table) && !is.data.frame(table)){
    stop("input x must be matrix or data frame\n")
  }

  pal <- colorRampPalette(c("white", col), space = "Lab")

  if(is.na(main)) {
    main <- deparse(substitute(table))
  }

  table <- as.matrix(table)

  nrow <- nrow(table)
  ncol <- ncol(table)

  m <- t(table)[, rev(seq(nrow))]

  image(m, main=main, xlab=xlab, ylab=ylab,
        xaxt=xaxt, yaxt=yaxt,
        col=pal(ncol*nrow), ...
  )

  grid(nx=ncol, ny=nrow, col=col)

  if(show.value) {
    text(rev(rep(seq(ncol), each=nrow)-1)/(ncol-1),
         (rep(seq(nrow), times=ncol)-1)/(nrow-1),
         rev(as.vector(table)), cex=2)
  }
}
