function Joint(x, y, numero) {
  this.xCentral = x;
  this.yCentral = y;
  this.x = x;
  this.y = y;
  this.size = 10;
  this.numero = str(numero);
  
  this.show = () => {
    strokeWeight(this.size)
    stroke(255);
    point(this.x, this.y);
    noStroke();
    fill(0);
    rectMode(CENTER);
  }
  
  this.update = (dx, dy) => {
    this.x = this.xCentral + dx;
    this.y = this.yCentral + dy; 
  }
}