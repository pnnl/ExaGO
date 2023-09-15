function LineColor(line) {
    var loading = Math.abs(line.properties.PF / line.properties.RATE_A);
    var r = Math.min(255, 255 * loading);
    var g = 0;
    var b = Math.max(0, 255 * (1 - loading));
  
    return [r, g, b];
  }
  
  function FlowColor(line){
    var loading = Math.abs(line.properties.PF / line.properties.RATE_A);
    var r = Math.min(255, 255 * loading);
    var g = 0;
    var b = Math.max(0, 255 * (1 - loading));
  
    return [r, g, b];
  }
  
  function FillColor(subst) {
    var Vm = subst.properties.bus[0].VM;
    var r = Math.min(255, 255 * (1.1 - Vm) / 0.2);
    var g = b;
    var b = Math.max(0, 255 * (Vm - 0.9) / 0.2);
  
    return [r, g, b];
  }
  
  
  function fillGenColumnColor(data) {
    if (data.color == 'red') return [255, 0, 0];
    else if (data.color == 'green') return [0, 255, 0];
    else if (data.color == 'yellow') return [244, 219, 135];
    else if (data.color == 'gray') return [128, 128, 128];
    else if (data.color == 'blue') return [28, 163, 236];
    else if (data.color == 'orange') return [255, 165, 0];
    else if (data.color == 'black') return [0, 0, 0];
  
  }
  
  function fillGenColumnColorCap(data) {
    var color;
  
    color = fillGenColumnColor(data);
  
    color = [...color, 255 * 0.3];
  
    return color;
  }
  
  
  function getVoltageFillColor(data) {
    var Vm = data.properties.Vm_avg;
    var r = Math.min(255, 255 * (1.1 - Vm) / 0.2);
    var b = 0;
    var g = Math.max(0, 255 * (Vm - 0.9) / 0.2);
  
    return [r, g, b];
  }

  export { LineColor, FlowColor, FillColor, fillGenColumnColor, fillGenColumnColorCap, getVoltageFillColor};
