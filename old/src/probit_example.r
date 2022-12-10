dev.new(width=10,height=7, units="in")
# png(filename="probit_example.png", width=10,height=7, units="in", res=300, type="cairo-png")
plot(dnorm(seq(-4,4,length=100)), seq(-4,4,length=100), type="l", lwd=3, col="darkgray", 
      xlab="", ylab="", xaxs="i", xaxt="n", yaxt="n", xlim=c(0,dnorm(0)*1.1*7))
mtext("Age", side=2, line=2, cex=1.5)
lines(c(0,dnorm(0)), c(0,0), col="blue", lwd=2)
arrows(dnorm(-1.5), -1.5, dnorm(1.5), 1.5, code=3, length=0.15, angle=20, col="blue", lwd=2)
# dev.off()

abline(h=c(-1.1,0.6), col="blue", lwd=2)
# dev.off()

mu.seq <- c(0,runif(6,-2,2))
for(i in 1:6)
{
  lines(dnorm(seq(-4,4,length=100), mean=mu.seq[i+1]) + i*dnorm(0)*1.1, 
        seq(-4,4,length=100), lwd=3, col="darkgray")
  abline(v=i*dnorm(0)*1.1)
}
mtext("Cohort", side=1, line=1, cex=1.5)
# dev.off()
  

