
## Toy model
############
trait_model <-
function(input=c(1,500,1,1,1,1)) {
  .C("trait_model",input=input,stat_to_return=array(0,4))$stat_to_return
}

