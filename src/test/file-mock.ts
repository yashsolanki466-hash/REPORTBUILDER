// Mocking File.text() since it's not fully implemented in some jsdom/node environments
if (!File.prototype.text) {
  File.prototype.text = async function () {
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = () => resolve(reader.result as string);
      reader.onerror = reject;
      reader.readAsText(this);
    });
  };
}
