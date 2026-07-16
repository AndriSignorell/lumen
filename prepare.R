

Rcpp::compileAttributes()
devtools::clean_dll()
devtools::document()
devtools::check()

.rs.restartR()
devtools::install()

devtools::build_manual(pkg = "C:/temp/DescToolsX")
devtools::build_manual(pkg = "C:/temp/lumen")
devtools::build_manual(pkg = "C:/temp/pharos")
devtools::build_manual(pkg = "C:/temp/bedrock")

devtools::document()
devtools::load_all()
devtools::check()
devtools::test()
devtools::run_examples()

covr::package_coverage()

available::available("figura", browse = FALSE)
available::available("playfair", browse = FALSE)
available::available("clarity", browse = FALSE)

pkgdown::build_site()


rename_rda <- function(old_file, old_name, new_name, data_dir = "C:/temp/bedrock/data/") {
  e <- new.env()
  load(file.path(data_dir, old_file), envir = e)
  assign(new_name, get(old_name, envir = e), envir = e)
  rm(list = old_name, envir = e)
  
  tmp <- tempfile(fileext = ".rda")
  save(list = new_name, envir = e, file = tmp, compress = "xz")
  
  file.remove(file.path(data_dir, old_file))
  file.copy(tmp, file.path(data_dir, paste0(new_name, ".rda")))
  file.remove(tmp)
  
  cat(sprintf("  %s -> %s\n", old_name, new_name))
}


rename_rda("d.pizza.rda",  "d.pizza",  "Pizza")
rename_rda("cards.rda",    "cards",    "Cards")
rename_rda("roulette.rda", "roulette", "Roulette")
rename_rda("tarot.rda",    "tarot",    "Tarot")


Prefix <- d.prefix
Units  <- d.units
usethis::use_data(Prefix, Units, internal = TRUE, overwrite = TRUE)


# put the datasets in sysdata.rda
usethis::use_data(d.prefix, d.units, internal = TRUE)



## barplot für helpfile grafik
{

  bplot <- function(i, pos) {
    b <- plotBar(unname(VADeaths[, i]), col="lightsteelblue", border=FALSE, yaxt="n", xlab="",
                 mar=c(1,1,4,0.5), grid=FALSE, main=gettextf("pos=%s", dQuote(pos)), 
                 cex.main=1, font.main=1, stamp=NA)
    barText(VADeaths[, i], b=b, horiz = FALSE, cex=1.1, border=NA, xpd=NA, pos = pos)
  }
  
  
  
  png(filename = "c:/temp/bartext.png", width=1000, height=1000/2, 
      res=120)
  {
    par(mfrow=c(1,4), xpd=NA)
    
    bplot(1, "topout")
    bplot(1, "topin")
    bplot(1, "mid")
    bplot(1, "bottomin")
  }
  dev.off()

}







x <- regPolygon(numVertices = 3)

convUnit

linScale

devtools::test()

colLegend()
example(colLegend)

plotFdist(bedrock::d.pizza$temperature, na.rm=TRUE)

DescToolsX::desc(bedrock::d.pizza$wrongpizza, na.rm=TRUE)


head(bedrock::d.pizza)



undebug(convUnit)
convUnit(100, "km/h", "m/s")

fm(0.00325, fmt="%")


plotViolin(len~ suppToothGrowth)

boxplot(len ~ supp : dose, ToothGrowth)

par(mfrow=c(2,1))
plotViolin(len ~ supp, ToothGrowth)

plotViolin(ToothGrowth$len, ylim=c(-10,50), box=TRUE, grid=list(col="red"))


plotViolin(len ~ supp:dose, ToothGrowth, grid=list(col="red"))



plotQQ(ToothGrowth$len, ylim=c(-10,50), grid=list(col="red"))
par("las")
bedrock::m

plotBinaryTree()

x <- sort(sample(100, 24))
z <- plotBinaryTree(x, cex=0.8)

drawRegPolygon()polygon

plotBinaryTree(LETTERS[1:15], 
               text=list(col="deeppink4", cex=2),
               line=list(col="navajowhite3", lwd=2))

devtools::build_manual(pkg = "C:/temp/aurora")
devtools::build_manual(pkg = "C:/temp/bedrock")
devtools::build_manual(pkg = "C:/temp/DescToolsX")
devtools::build_manual(pkg = "C:/temp/lumen")


mtcars
mtcars

normalizePath(file.path(Sys.getenv("APPDATA"),
                        "RStudio",
                        "keybindings"))

# C:\\Users\\Andri\\AppData\\Roaming\\RStudio\\keybindings





plotBar(VADeaths/1e3,  box=FALSE, bg="lightyellow", main="VADeaths",
        horiz=TRUE, 
        text=list(border=FALSE, cex=0.8, col=c("blue", "green","orange")), 
        mar=c(right=5), 
        yax = list(fmt="%", d=0, big=",", 
                   col="red", col.axis="blue", lwd=2))



debug(pal)
pal()

preview()
